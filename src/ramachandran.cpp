#include "ramachandran.h"

#include <core/log.h>
#include <core/math_utils.h>

#include <mol/molecule_dynamic.h>
#include <mol/molecule_utils.h>
#include <mol/trajectory_utils.h>

#include "gfx/gl.h"
#include "gfx/gl_utils.h"
#include "image.h"
#include "task_system.h"

struct Coord {
    i16 x, y;
};

struct BackboneAnglesTrajectory {
    i32 num_segments = 0;
    i32 num_frames = 0;
    Array<BackboneAngle> angle_data{};
};

static inline Array<BackboneAngle> get_frame_backbone_angles(BackboneAnglesTrajectory& backbone_angle_traj, int frame_index) {
    if (backbone_angle_traj.angle_data.count == 0 || backbone_angle_traj.num_segments == 0) return {};
    ASSERT(frame_index < backbone_angle_traj.angle_data.count / backbone_angle_traj.num_segments);
    return Array<BackboneAngle>(&backbone_angle_traj.angle_data[frame_index * backbone_angle_traj.num_segments], backbone_angle_traj.num_segments);
}

static inline Array<BackboneAngle> get_range_backbone_angles(BackboneAnglesTrajectory& backbone_angle_traj, Range<i32> range) {
    if (backbone_angle_traj.angle_data.count == 0 || backbone_angle_traj.num_segments == 0) return {};
    ASSERT(range.beg < backbone_angle_traj.num_frames);
    ASSERT(range.end <= backbone_angle_traj.num_frames);
    return backbone_angle_traj.angle_data.subarray(range.beg * backbone_angle_traj.num_segments, range.ext() * backbone_angle_traj.num_segments);
}

static inline i32 get_backbone_angles_trajectory_current_frame_count(const BackboneAnglesTrajectory& backbone_angle_traj) {
    if (backbone_angle_traj.angle_data.count == 0 || backbone_angle_traj.num_segments == 0) return 0;
    return (i32)backbone_angle_traj.angle_data.count / backbone_angle_traj.num_segments;
}

static inline Array<BackboneAngle> get_backbone_angles(BackboneAnglesTrajectory& backbone_angle_traj, int frame_index, AtomRange range) {
    return get_frame_backbone_angles(backbone_angle_traj, frame_index).subarray(range);
}

static void init_backbone_angles_trajectory(BackboneAnglesTrajectory* data, const MoleculeDynamic& dynamic);
static void free_backbone_angles_trajectory(BackboneAnglesTrajectory* data);
static task_system::ID compute_backbone_angles_trajectory(BackboneAnglesTrajectory* data, const MoleculeDynamic& dynamic);

namespace ramachandran {

static BackboneAnglesTrajectory traj_angles = {};

// Accumulation texture data
constexpr int acc_width = 1024;
constexpr int acc_height = 1024;

static GLuint gui_tex = 0;
static GLuint seg_tex = 0;
static GLuint acc_tex = 0;
static GLuint col_tex = 0;

static GLuint coord_buf = 0;
static GLuint fbo = 0;
static GLuint program = 0;

static GLuint vao = 0;
static GLuint vbo = 0;

static GLint uniform_loc_radius = -1;
static GLint uniform_loc_color = -1;

static Image src_img = {};  // This is the unmodified source image (Ramachandran plot)
static Image gui_img = {};  // Visual representation which is used in gui-elements
static Image seg_img = {};  // Segmentation version of image (Blurred, to hide low res artifacts and for smoother transitions)
static Image col_img = {};  // Color version of image (Even more blurred, for smoother transitions of secondary structure colors)

GLuint get_accumulation_texture() { return acc_tex; }
GLuint get_gui_texture() { return gui_tex; }
GLuint get_segmentation_texture() { return seg_tex; }
GLuint get_color_texture() { return col_tex; }

const Image& get_gui_image() { return gui_img; }
const Image& get_segmentation_image() { return seg_img; }
const Image& get_color_image() { return col_img; }

// @NOTE: This should generate a quad with a certain size in texture coordinates
constexpr const char* v_shader_src = R"(
#version 150 core

uniform float u_radius;

in vec2 pos;

void main() {
	gl_PointSize = u_radius * 2.0;
	gl_Position = vec4(pos, 0, 1);
}
)";

// @NOTE: Do some radial falloff based on uv coordinate
constexpr const char* f_shader_src = R"(
#version 150 core

uniform vec4 u_color;

out vec4 out_frag;

void main() {
    vec2 uv = gl_PointCoord.xy * 2.0 - 1.0;
	float dist2 = dot(uv, uv);
	float falloff = max(0.0, 1.0 - dist2);
    out_frag = vec4(u_color.rgb, u_color.a * falloff);
}
)";

void init_map(Image* img, GLuint tex, const ColorMap& color_map, int blur_level) {
    ASSERT(src_img);
    ASSERT(img);
    ASSERT(tex);

    init_image(img, src_img);

    constexpr u32 IN_BACKGROUND = 0xFFFFFFFF;
    constexpr u32 IN_ALPHA_HIGH = 0xFF0000FF;
    constexpr u32 IN_ALPHA_MID = 0xFF7F7FFF;
    constexpr u32 IN_BETA_HIGH = 0xFFFF0000;
    constexpr u32 IN_BETA_MID = 0xFFFF7F7F;
    constexpr u32 IN_LEFT_ALPHA_HIGH = 0xFF00FF00;
    constexpr u32 IN_LEFT_ALPHA_MID = 0xFF7FFF7F;
    constexpr u32 IN_P_MID = 0xFF7FFFFF;

    for (int i = 0; i < src_img.width * src_img.height; i++) {
        const u32 in_color = src_img.data[i];
        u32& out_color = img->data[i];
        switch (in_color) {
            case IN_BACKGROUND:
                out_color = math::convert_color(color_map.region_color[Region_None]);
                break;
            case IN_ALPHA_HIGH:
                out_color = math::convert_color(color_map.region_color[Region_AlphaHigh]);
                break;
            case IN_ALPHA_MID:
                out_color = math::convert_color(color_map.region_color[Region_AlphaMid]);
                break;
            case IN_BETA_HIGH:
                out_color = math::convert_color(color_map.region_color[Region_BetaHigh]);
                break;
            case IN_BETA_MID:
                out_color = math::convert_color(color_map.region_color[Region_BetaMid]);
                break;
            case IN_LEFT_ALPHA_HIGH:
                out_color = math::convert_color(color_map.region_color[Region_LeftAlphaHigh]);
                break;
            case IN_LEFT_ALPHA_MID:
                out_color = math::convert_color(color_map.region_color[Region_LeftAlphaMid]);
                break;
            case IN_P_MID:
                out_color = math::convert_color(color_map.region_color[Region_PMid]);
                break;
            default:
                break;
        }
    }

    gaussian_blur(img, blur_level);

    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, img->width, img->height, 0, GL_RGBA, GL_UNSIGNED_BYTE, img->data);
    glBindTexture(GL_TEXTURE_2D, 0);
}

void init_gui_map(const ColorMap& color_map, int blur_level) { init_map(&gui_img, gui_tex, color_map, blur_level); }
void init_segmentation_map(const ColorMap& color_map, int blur_level) { init_map(&seg_img, seg_tex, color_map, blur_level); }
void init_color_map(const ColorMap& color_map, int blur_level) { init_map(&col_img, col_tex, color_map, blur_level); }

static void update_vbo();

task_system::ID initialize(const MoleculeDynamic& dynamic) {
    task_system::ID id = task_system::INVALID_ID;

    if (dynamic) {
        init_backbone_angles_trajectory(&traj_angles, dynamic);
        id = compute_backbone_angles_trajectory(&traj_angles, dynamic);
        task_system::wait_for_task(id);
        update_vbo();
    }

    if (!read_image(&src_img, VIAMD_IMAGE_DIR "/ramachandran.bmp")) {
        LOG_ERROR("Could not read ramachandran map!");
        return id;
    }

    if (!gui_tex) {
        glGenTextures(1, &gui_tex);
        glBindTexture(GL_TEXTURE_2D, gui_tex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    if (!seg_tex) {
        glGenTextures(1, &seg_tex);
        glBindTexture(GL_TEXTURE_2D, seg_tex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    if (!col_tex) {
        glGenTextures(1, &col_tex);
        glBindTexture(GL_TEXTURE_2D, col_tex);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    {
        ColorMap seg_color_map;
        seg_color_map.region_color[Region_None] = vec4(0, 0, 0, 0);
        seg_color_map.region_color[Region_AlphaHigh] = vec4(1, 0, 0, 1);
        seg_color_map.region_color[Region_AlphaMid] = vec4(0, 0, 0, 0);
        seg_color_map.region_color[Region_BetaHigh] = vec4(0, 0, 1, 1);
        seg_color_map.region_color[Region_BetaMid] = vec4(0, 0, 0, 0);
        seg_color_map.region_color[Region_LeftAlphaHigh] = vec4(1, 0, 0, 1);
        seg_color_map.region_color[Region_LeftAlphaMid] = vec4(0, 0, 0, 0);
        seg_color_map.region_color[Region_PMid] = vec4(0, 0, 0, 0);
        init_segmentation_map(seg_color_map, 2);
    }

    constexpr float h_r = 0.0f / 360.0f;
    constexpr float h_g = 120.0f / 360.0f;
    constexpr float h_b = 240.0f / 360.0f;
    constexpr float h_y = 60.0f / 360.0f;

    {
        constexpr float c = 0.45f;
        constexpr float l = 0.9f;
        constexpr float a_high = 0.8f;
        constexpr float a_mid = 0.2f;

        ColorMap gui_color_map;
        gui_color_map.region_color[Region_None] = vec4(0, 0, 0, 0);
        gui_color_map.region_color[Region_AlphaHigh] = vec4(math::hcl_to_rgb(h_r, c, l), a_high);
        gui_color_map.region_color[Region_AlphaMid] = vec4(math::hcl_to_rgb(h_r, c, l), a_mid);
        gui_color_map.region_color[Region_BetaHigh] = vec4(math::hcl_to_rgb(h_b, c, l), a_high);
        gui_color_map.region_color[Region_BetaMid] = vec4(math::hcl_to_rgb(h_b, c, l), a_mid);
        gui_color_map.region_color[Region_LeftAlphaHigh] = vec4(math::hcl_to_rgb(h_g, c, l), a_high);
        gui_color_map.region_color[Region_LeftAlphaMid] = vec4(math::hcl_to_rgb(h_g, c, l), a_mid);
        gui_color_map.region_color[Region_PMid] = vec4(math::hcl_to_rgb(h_y, c, l), a_mid);
        init_gui_map(gui_color_map, 2);
    }

    {
        constexpr float c = 0.8f;
        constexpr float l = 0.9f;
        constexpr float a_high = 0.8f;
        constexpr float a_mid = 0.2f;

        ColorMap col_color_map;
        col_color_map.region_color[Region_None] = vec4(1, 1, 1, 1);
        col_color_map.region_color[Region_AlphaHigh] = vec4(math::hcl_to_rgb(h_r, c, l), a_high);
        col_color_map.region_color[Region_AlphaMid] = vec4(math::hcl_to_rgb(h_r, c, l), a_mid);
        col_color_map.region_color[Region_BetaHigh] = vec4(math::hcl_to_rgb(h_b, c, l), a_high);
        col_color_map.region_color[Region_BetaMid] = vec4(math::hcl_to_rgb(h_b, c, l), a_mid);
        col_color_map.region_color[Region_LeftAlphaHigh] = vec4(math::hcl_to_rgb(h_g, c, l), a_high);
        col_color_map.region_color[Region_LeftAlphaMid] = vec4(math::hcl_to_rgb(h_g, c, l), a_mid);
        col_color_map.region_color[Region_PMid] = vec4(math::hcl_to_rgb(h_y, c, l), a_mid);
        init_color_map(col_color_map, 4);
    }

    if (!program) {
        GLuint v_shader = gl::compile_shader_from_source(v_shader_src, GL_VERTEX_SHADER);
        GLuint f_shader = gl::compile_shader_from_source(f_shader_src, GL_FRAGMENT_SHADER);
        defer {
            glDeleteShader(v_shader);
            glDeleteShader(f_shader);
        };

        program = glCreateProgram();
        const GLuint shaders[] = {v_shader, f_shader};
        gl::attach_link_detach(program, shaders);

        uniform_loc_radius = glGetUniformLocation(program, "u_radius");
        uniform_loc_color = glGetUniformLocation(program, "u_color");
    }


    if (!acc_tex) {
        glGenTextures(1, &acc_tex);
        glBindTexture(GL_TEXTURE_2D, acc_tex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, acc_width, acc_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);
    }

    if (!coord_buf) {
        glGenBuffers(1, &coord_buf);
    }

    if (!fbo) {
        glGenFramebuffers(1, &fbo);
        glBindFramebuffer(GL_FRAMEBUFFER, fbo);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, acc_tex, 0);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    if (!vao) {
        glGenVertexArrays(1, &vao);
    }

    return id;
}

void shutdown() {
    free_backbone_angles_trajectory(&traj_angles);

    if (seg_tex) glDeleteTextures(1, &seg_tex);
    if (acc_tex) glDeleteTextures(1, &acc_tex);
    if (col_tex) glDeleteTextures(1, &col_tex);
    if (coord_buf) glDeleteBuffers(1, &coord_buf);
    if (fbo) glDeleteFramebuffers(1, &fbo);
    if (vbo) glDeleteBuffers(1, &vbo);
    if (vbo) glDeleteVertexArrays(1, &vao);
}

void clear_accumulation_texture() {
    GLint last_viewport[4];
    glGetIntegerv(GL_VIEWPORT, last_viewport);

    glViewport(0, 0, acc_width, acc_height);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);

    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT);

    glViewport(last_viewport[0], last_viewport[1], (GLsizei)last_viewport[2], (GLsizei)last_viewport[3]);
}

void render_accumulation_texture(Range<i32> frame_range, vec4 color, float radius) {
    // Backup GL state
    GLint last_polygon_mode[2];
    glGetIntegerv(GL_POLYGON_MODE, last_polygon_mode);
    GLint last_viewport[4];
    glGetIntegerv(GL_VIEWPORT, last_viewport);
    GLenum last_blend_src_rgb;
    glGetIntegerv(GL_BLEND_SRC_RGB, (GLint*)&last_blend_src_rgb);
    GLenum last_blend_dst_rgb;
    glGetIntegerv(GL_BLEND_DST_RGB, (GLint*)&last_blend_dst_rgb);
    GLenum last_blend_src_alpha;
    glGetIntegerv(GL_BLEND_SRC_ALPHA, (GLint*)&last_blend_src_alpha);
    GLenum last_blend_dst_alpha;
    glGetIntegerv(GL_BLEND_DST_ALPHA, (GLint*)&last_blend_dst_alpha);
    GLenum last_blend_equation_rgb;
    glGetIntegerv(GL_BLEND_EQUATION_RGB, (GLint*)&last_blend_equation_rgb);
    GLenum last_blend_equation_alpha;
    glGetIntegerv(GL_BLEND_EQUATION_ALPHA, (GLint*)&last_blend_equation_alpha);
    GLboolean last_enable_blend = glIsEnabled(GL_BLEND);
    GLboolean last_enable_cull_face = glIsEnabled(GL_CULL_FACE);
    GLboolean last_enable_depth_test = glIsEnabled(GL_DEPTH_TEST);
    GLboolean last_enable_scissor_test = glIsEnabled(GL_SCISSOR_TEST);

    const GLint offset = frame_range.beg * traj_angles.num_segments;
    const GLsizei count = (GLsizei)frame_range.ext() * traj_angles.num_segments;

    // RENDER TO ACCUMULATION TEXTURE
    glViewport(0, 0, acc_width, acc_height);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, fbo);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);

    glEnable(GL_BLEND);
    glBlendEquation(GL_FUNC_ADD);
    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_PROGRAM_POINT_SIZE);

    glUseProgram(program);

    const float rad = math::max(radius, 0.5f);
    const vec4  col = color * vec4(1.0f, 1.0f, 1.0f, 1.0f - (rad - radius));

    // Draw
    glUniform1f(uniform_loc_radius, rad);
    glUniform4fv(uniform_loc_color, 1, &col[0]);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glVertexAttribPointer(0, 2, GL_SHORT, GL_TRUE, sizeof(Coord), (const GLvoid*)0);
    glEnableVertexAttribArray(0);

    glDrawArrays(GL_POINTS, offset, count);

    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    glUseProgram(0);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);

    glDisable(GL_PROGRAM_POINT_SIZE);

    // Restore modified GL state
    glBlendEquationSeparate(last_blend_equation_rgb, last_blend_equation_alpha);
    glBlendFuncSeparate(last_blend_src_rgb, last_blend_dst_rgb, last_blend_src_alpha, last_blend_dst_alpha);
    if (last_enable_blend)
        glEnable(GL_BLEND);
    else
        glDisable(GL_BLEND);
    if (last_enable_cull_face)
        glEnable(GL_CULL_FACE);
    else
        glDisable(GL_CULL_FACE);
    if (last_enable_depth_test)
        glEnable(GL_DEPTH_TEST);
    else
        glDisable(GL_DEPTH_TEST);
    if (last_enable_scissor_test)
        glEnable(GL_SCISSOR_TEST);
    else
        glDisable(GL_SCISSOR_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, last_polygon_mode[0]);
    glViewport(last_viewport[0], last_viewport[1], (GLsizei)last_viewport[2], (GLsizei)last_viewport[3]);
}

void update_vbo() {
    if (traj_angles.angle_data.size() == 0) return;

    if (!vbo) {
        glGenBuffers(1, &vbo);
    }

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, traj_angles.angle_data.size() * sizeof(Coord), NULL, GL_STATIC_DRAW);

    void* ptr = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    if (!ptr) {
        LOG_ERROR("Could not map angle buffer");
        return;
    }

    Coord* coords = (Coord*)ptr;
    constexpr float ONE_OVER_PI = 1.f / math::PI;
    for (i64 i = 0; i < traj_angles.angle_data.size(); i++) {
        const BackboneAngle& angle = traj_angles.angle_data[i];
        if (angle.phi == 0 || angle.psi == 0) continue;
        vec2 coord = vec2(angle.phi, -angle.psi) * ONE_OVER_PI;  // [-PI, PI] -> [-1, 1]
        //coord = (vec2(ivec2(coord * vec2(acc_width / 2, acc_height / 2)))) / vec2(acc_width / 2, acc_height / 2);
        coords[i].x = (short)(coord.x * INT16_MAX); 
        coords[i].y = (short)(coord.y * INT16_MAX);
    }

    glUnmapBuffer(GL_ARRAY_BUFFER);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

}  // namespace ramachandran

static void init_backbone_angles_trajectory(BackboneAnglesTrajectory* data, const MoleculeDynamic& dynamic) {
    ASSERT(data);
    if (!dynamic.molecule || !dynamic.molecule.residue.backbone.atoms || !dynamic.trajectory) return;

    if (data->angle_data) {
        FREE(data->angle_data.ptr);
    }

    i32 alloc_count = (i32)dynamic.molecule.residue.count * (i32)dynamic.trajectory.frame_buffer.count;
    data->num_segments = (i32)dynamic.molecule.chain.count;
    data->num_frames = 0;
    data->angle_data = {(BackboneAngle*)CALLOC(alloc_count, sizeof(BackboneAngle)), alloc_count};
}

static void free_backbone_angles_trajectory(BackboneAnglesTrajectory* data) {
    ASSERT(data);
    if (data->angle_data) {
        FREE(data->angle_data.ptr);
        *data = {};
    }
}

static task_system::ID compute_backbone_angles_trajectory(BackboneAnglesTrajectory* data, const MoleculeDynamic& dynamic) {
    ASSERT(dynamic);
    if (dynamic.trajectory.num_frames == 0 || !dynamic.molecule.residue.backbone.angle) return task_system::INVALID_ID;

    task_system::ID compute_task = task_system::enqueue_pool(
        "Backbone Angles Trajectory", dynamic.trajectory.num_frames, [data, &dynamic](task_system::TaskSetRange range) {
            for (u32 f_idx = range.beg; f_idx < range.end; f_idx++) {
                const soa_vec3 pos = get_trajectory_positions(dynamic.trajectory, f_idx);
                Array<BackboneAngle> frame_angles = get_frame_backbone_angles(*data, f_idx);
                for (i64 ci = 0; ci < dynamic.molecule.chain.count; ++ci) {
                    const auto bb_range = dynamic.molecule.chain.residue_range[ci];
                    const auto bb_seg = dynamic.molecule.residue.backbone.atoms + bb_range.beg;
                    auto bb_ang = frame_angles.data() + bb_range.beg;
                    if (bb_range.ext() < 2) {
                        *bb_ang = {};
                    } else {
                        compute_backbone_angles(bb_ang, pos, bb_seg, bb_range.ext());
                    }
                }
            }
        });

    return compute_task;
}