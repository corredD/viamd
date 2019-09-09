#include "property.h"

#include <core/math_utils.h>
#include <core/hash.h>
#include <core/log.h>
#include <mol/molecule_structure.h>
#include <mol/molecule_trajectory.h>
#include <mol/molecule_dynamic.h>
#include <mol/molecule_utils.h>
#include <mol/trajectory_utils.h>

#include <thread>
#include <tinyexpr.h>

#define COMPUTE_ID(x) (hash::crc64(x))

namespace prop {

typedef bool (*StructureMatchFunc)(Bitfield out_atoms, String out_err_msg, Array<CString> args, const MoleculeStructure& mol);
typedef bool (*PropertyComputeFunc)(Array<float> out_data, String out_err_msg, Array<CString> args, const MoleculeTrajectory& traj);

struct PropertyComputeCommand {
    ID id = INVALID_ID;
    PropertyComputeFunc func = nullptr;
};

struct StructureMatchEntry {
    ID id = INVALID_ID;
    StructureMatchFunc func = nullptr;
};

struct Property {
    ID id;
    Array<float> data;
    bool periodic;
    bool valid;

    StringBuffer<32> name;
    StringBuffer<32> unit;
    StringBuffer<1024> error_msg;
};

struct StatisticsContext {
    DynamicArray<PropertyComputeCommand> property_func_entries{};
    DynamicArray<StructureMatchEntry> structure_func_entries{};
    DynamicArray<Property*> properties{};

    Property* current_property = nullptr;

    volatile bool thread_running = false;
    volatile bool stop_signal = false;
    volatile float fraction_done = 0.f;
};

static StatisticsContext ctx;

static PropertyComputeFunc find_property_compute_func(CString cmd) {
    const ID id = COMPUTE_ID(cmd);
    for (auto& e : ctx.property_func_entries) {
        if (id == e.id) {
            return e.func;
        }
    }
    return nullptr;
}

static StructureMatchFunc find_structure_match_func(CString cmd) {
    const ID id = COMPUTE_ID(cmd);
    for (auto& e : ctx.structure_func_entries) {
        if (id == e.id) {
            return e.func;
        }
    }
    return nullptr;
}

static Range<float> compute_range(Array<const float> data) {
    if (data.count == 0) {
        return {0, 0};
    }
    Range<float> range{FLT_MAX, -FLT_MAX};
    for (float v : data) {
        range.beg = math::min(range.beg, v);
        range.end = math::max(range.end, v);
    }
    return range;
}

static DynamicArray<CString> extract_arguments(CString str) {
    DynamicArray<CString> args;

    const auto* beg = str.beg();
    const auto* end = str.beg();
    int32 count = 0;

    while (end < str.end()) {
        if (*end == '(')
            count++;
        else if (*end == ')')
            count--;
        else if (*end == ',') {
            if (count == 0) {
                args.push_back(trim({beg, end}));
                beg = end + 1;
            }
        }
        end++;
    }
    if (beg != end) args.push_back(trim({beg, end}));

    return args;
}

static CString extract_command(CString str) {
    str = trim(str);
    const char* ptr = str.beg();
    while (ptr != str.end() && *ptr != '(' && !is_whitespace(*ptr)) ptr++;
    return {str.beg(), ptr};
}

void set_error_message(String dst, const char* format, ...) {
    va_list ap;
    va_start(ap, format);
    int count = vsnprintf(dst.data(), dst.size(), format, ap);
    LOG_ERROR(format, ap);
    va_end(ap);
}

bool structure_match_residue(Bitfield atoms, String err_msg, const Array<CString> args, const MoleculeStructure& mol) {
    // Expect args to be  > 0
    if (args.count == 0) {
        snprintf(err_msg.data(), err_msg.size(), "Expects one or more arguments for residue");
        return false;
    }

    for (const auto& arg : args) {
        Range<int32> range = {-1, -1};

        if (is_alpha(arg[0])) {
        } else if (is_range(arg)) {
            if (!extract_range(&range, arg)) {
                set_error_message(err_msg, "Failed to parse range in argument for residue");
                return false;
            }
            if (range.x == -1) range.x = 1;
            if (range.y == -1) range.y = (int32)mol.residues.count;
        } else if (const auto id = to_int(arg)) {
            if (!id.success) {
                set_error_message(err_msg, "Failed to parse argument for residue");
                return false;
            }
            range.x = range.y = id;
        } else {
            set_error_message(err_msg, "Could not parse argument in residue: '%.s'", arg.length(), arg.cstr());
        }

        if (range.x < 1 || (int32)mol.residues.count < range.y) {
            set_error_message(err_msg, "Index for residue is out of bounds");
            return false;
        }
        for (int32 i = range.x - 1; i < range.y; i++) {
            const auto& res = mol.residues[i];
            bitfield::set_range(atoms, res.atom_range);
        }
    }

    return true;
}

bool structure_match_chain(Bitfield atoms, String err_msg, const Array<CString> args, const MoleculeStructure& mol) {
    // Expect args.count to be > 0
    if (args.count == 0) {
        set_error_message(err_msg, "Expects one or more arguments for chain");
        return false;
    }

    for (const auto& chain : mol.chains) {
        for (const auto& arg : args) {
            if (compare(chain.id, arg)) {
                bitfield::set_range(atoms, chain.atom_range);
                break;
            }
        }
    }

    return true;
}

#include "rmsd.h"

static float rmsd(const float* ref_x, const float* ref_y, const float* ref_z, const float* cur_x, const float* cur_y, const float* cur_z, int64 count) {
    if (count <= 1) return 0.f;

    // ugly ugly hacks
    double* ref_tmp = (double*)TMP_MALLOC(count * sizeof(double) * 3);
    double* cur_tmp = (double*)TMP_MALLOC(count * sizeof(double) * 3);
    defer {
        TMP_FREE(ref_tmp);
        TMP_FREE(cur_tmp);
    };
    for (int64 i = 0; i < count; i++) {
        ref_tmp[i * 3 + 0] = ref_x[i];
        ref_tmp[i * 3 + 1] = ref_y[i];
        ref_tmp[i * 3 + 2] = ref_z[i];
        cur_tmp[i * 3 + 0] = cur_x[i];
        cur_tmp[i * 3 + 1] = cur_y[i];
        cur_tmp[i * 3 + 2] = cur_z[i];
    }

    double val;
    fast_rmsd((double(*)[3])ref_tmp, (double(*)[3])cur_tmp, (int)count, &val);

    return (float)val;
}

static DynamicArray<Property*> extract_property_dependencies(Array<Property*> properties, CString expression) {
    DynamicArray<Property*> dependencies;
    for (Property* prop : properties) {
        if (!prop->valid) continue;
        CString match = find_string(expression, prop->name);
        if (match) {
            if (match.beg() != expression.beg()) {
                char pre_char = *(match.beg() - 1);
                if (is_alpha(pre_char)) continue;
                if (is_digit(pre_char)) continue;
            }
            if (match.end() != expression.end()) {
                char post_char = *match.end();
                if (is_alpha(post_char)) continue;
                if (is_digit(post_char)) continue;
            }
            dependencies.push_back(prop);
        }
    }
    return dependencies;
}

/*
static bool compute_expression(Property* prop, const Array<CString> args, const MoleculeDynamic&) {
    ASSERT(prop);
    if (args.count == 0) {
        set_error_message(err_msg, "expression expects 1 or more arguments");
        return false;
    }

    // Concatenate all arguments
    StringBuffer<1024> expr_str = CString(args.front().beg(), args.back().end());

    if (!expr_str) {
        return false;
    }

    if (!balanced_parentheses(expr_str)) {
        set_error_message(err_msg, "Expression contains unbalanced parentheses!");
        return false;
    }

    // Extract which properties preceedes this property
    Array<Property*> properties = ctx.properties;
    properties.count = 0;
    for (int i = 0; i < ctx.properties.size(); i++) {
        if (prop == ctx.properties[i]) {
            properties.count = i;
            break;
        }
    }

    prop->dependencies = extract_property_dependencies(properties, expr_str);

    DynamicArray<double> values(prop->dependencies.size(), 0);
    DynamicArray<te_variable> vars;
    for (int32 i = 0; i < prop->dependencies.size(); i++) {
        vars.push_back({prop->dependencies[i]->name_buf.cstr(), &values[i], 0, 0});
    }

    int err;
    te_expr* expr = te_compile(expr_str.cstr(), vars.data(), (int32)vars.size(), &err);

    if (expr) {
        int32 max_instance_count = 0;
        for (auto* p : prop->dependencies) {
            max_instance_count = math::max(max_instance_count, (int32)p->instance_data.count);
        }

        const int32 frame_count = (int32)prop->avg_data.size();
        init_instance_data(&prop->instance_data, max_instance_count, frame_count);

        float scl = 1.f / (float)max_instance_count;
        for (int32 frame = 0; frame < frame_count; frame++) {
            float val = 0.f;
            for (int32 i = 0; i < max_instance_count; i++) {
                for (int32 j = 0; j < values.size(); j++) {
                    values[j] = 0;
                    if (prop->dependencies[j]->instance_data.count == max_instance_count) {
                        if (frame < prop->dependencies[j]->instance_data[i].data.count) {
                            values[j] = prop->dependencies[j]->instance_data[i].data[frame];
                        }
                    } else {
                        if (frame < prop->dependencies[j]->instance_data[0].data.count) {
                            values[j] = prop->dependencies[j]->instance_data[0].data[frame];
                        }
                    }
                }
                prop->instance_data[i].data[frame] = (float)te_eval(expr);
                val += prop->instance_data[i].data[frame];
            }
            prop->avg_data[frame] = val * scl;
        }

        te_free(expr);
    } else {
        set_error_message(err_msg, "Malformed expression:\n%s\n%*s^\nError near here", expr_str.cstr(), err - 1, "");
        return false;
    }

    prop->total_data_range = compute_range(*prop);
    prop->avg_data_range = compute_range(prop->avg_data);
    prop->periodic = false;
    prop->unit_buf = "";

    return true;
}
*/

bool extract_atoms_from_argument(Bitfield atom_field, bool& com, CString arg) {
    
}

vec3 extract_mean_position_from_field(const float* x, const float* y, const float* z, const Bitfield field) {
    vec3 pos = {0, 0, 0};
    int64 i = bitfield::find_next_bit_set(field);
    int64 count = 0;
    while (i != -1) {
        pos.x += x[i];
        pos.y += y[i];
        pos.z += z[i];
    };
    i = bitfield::find_next_bit_set(field, i + 1);
    ++count;

    if (count > 0) {
        pos /= (float)count;
    }

    return pos;
}  // namespace prop

void extract_positions_from_field(Array<vec3> dst_pos, const float* x, const float* y, const float* z, const Bitfield field) {
    int64 src_idx = bitfield::find_next_bit_set(field);
    int64 dst_idx = 0;

    while (src_idx != -1 && dst_idx < dst_pos.size()) {
        dst_pos[dst_idx] = {x[src_idx], y[src_idx], z[src_idx]};
        src_idx = bitfield::find_next_bit_set(field, src_idx + 1);
    }
}

bool compute_distance(Array<float> out_data, String err_buf, Array<CString> args, const MoleculeTrajectory& traj) {
    if (args.size() != 2) {
        set_error_message(err_buf, "distance expects exactly two arguments");
        return false;
    }

    Bitfield field[2];
    bitfield::init(&field[0], out_data.size());
    bitfield::init(&field[1], out_data.size());
    defer {
        bitfield::free(&field[0]);
        bitfield::free(&field[1]);
    };

    bool com[2];
    if (!extract_atoms_from_argument(field[0], com[0], args[0])) {
        set_error_message(err_buf, "Could not extract atomic indices from first argument");
        return false;
    }

    if (!extract_atoms_from_argument(field[1], com[1], args[1])) {
        set_error_message(err_buf, "Could not extract atomic indices from second argument");
        return false;
    }

    vec3 pos[2];
    for (int64 i = 0; i < out_data.size(); ++i) {
        const auto x = get_trajectory_position_x(traj, i);
        const auto y = get_trajectory_position_y(traj, i);
        const auto z = get_trajectory_position_z(traj, i);

        if (com[0]) {
            pos[0] = extract_mean_position_from_field(x.data(), y.data(), z.data(), field[0]);
        } else {
            extract_positions_from_field({&pos[0], 1}, x.data(), y.data(), z.data(), field[0]);
        }

        if (com[1]) {
            pos[1] = extract_mean_position_from_field(x.data(), y.data(), z.data(), field[1]);
        } else {
            extract_positions_from_field({&pos[1], 1}, x.data(), y.data(), z.data(), field[1]);
        }

        out_data[i] = math::distance(pos[0], pos[1]);
    }

    return true;
}

void initialize() {
    ctx.property_func_entries.push_back({COMPUTE_ID("distance"), compute_distance});
    ctx.property_func_entries.push_back({COMPUTE_ID("angle"), compute_angle});
    ctx.property_func_entries.push_back({COMPUTE_ID("dihedral"), compute_dihedral});
    ctx.property_func_entries.push_back({COMPUTE_ID("rmsd"), compute_rmsd});
    ctx.property_func_entries.push_back({COMPUTE_ID("expression"), compute_expression});

    ctx.structure_func_entries.push_back({COMPUTE_ID("resname"), structure_match_resname});
    ctx.structure_func_entries.push_back({COMPUTE_ID("resid"), structure_match_resid});
    ctx.structure_func_entries.push_back({COMPUTE_ID("residue"), structure_match_residue});
    ctx.structure_func_entries.push_back({COMPUTE_ID("chainid"), structure_match_chainid});
    ctx.structure_func_entries.push_back({COMPUTE_ID("chain"), structure_match_chain});
    ctx.structure_func_entries.push_back({COMPUTE_ID("atom"), structure_match_atom});
    ctx.structure_func_entries.push_back({COMPUTE_ID("resatom"), structure_extract_resatom});
    ctx.structure_func_entries.push_back({COMPUTE_ID("com"), structure_apply_aggregation_strategy_com});
}

void shutdown() {}

bool register_property_command(CString cmd_keyword, PropertyComputeFunc compute_func, PropertyVisualizeFunc visualize_func) {
    if (!cmd_keyword.ptr || cmd_keyword.count == 0) {
        LOG_ERROR("Property command cannot be an empty string!");
        return false;
    }

    if (contains_whitespace(cmd_keyword)) {
        LOG_ERROR("Property command cannot contain whitespace!");
        return false;
    }

    if (!compute_func) {
        LOG_ERROR("Property command must have compute function!");
        return false;
    }

    ID hash = COMPUTE_ID(cmd_keyword);
    if (find_property_func_entry(hash) != nullptr) {
        LOG_ERROR("Property command already registered!");
        return false;
    }

    ctx.property_func_entries.push_back({hash, compute_func, visualize_func});
    return true;
}

bool sync_structure_data_length(Array<StructureData> data) {
    int32 max_count = 0;
    for (const auto& s : data) {
        max_count = math::max(max_count, (int32)s.structures.size());
    }

    // Extend and copy data from first element to match multi_count
    for (auto& s : data) {
        while (s.structures.size() < max_count) {
            s.structures.push_back(s.structures.front());
        }
    }

    if (data.count > 0) {
        for (int32 i = 0; i < data[0].structures.size(); i++) {
            int32 max_structure_count = 0;
            for (const auto& s : data) {
                const int32 c = s.strategy == COM ? 1 : structure_index_count(s.structures[i]);
                max_structure_count = math::max(max_structure_count, c);
                if (c > 1 && c != max_structure_count) {
                    set_error_message(err_msg, "Structures matched has different sizes in different arguments, this is not supported");
                    return false;
                }
            }
        }
    }

    return true;
}

static bool compute_property_data(Property* prop, const MoleculeDynamic& dynamic, int32 num_frames) {
    ASSERT(prop);

    ctx.current_property = prop;
    prop->error_msg_buf = "";
    if (prop->avg_data.size() != num_frames) {
        prop->avg_data.resize(num_frames);
    }
    if (prop->std_dev_data.size() != num_frames) {
        prop->std_dev_data.resize(num_frames);
    }
    if (prop->filter_fraction.size() != num_frames) {
        prop->filter_fraction.resize(num_frames);
    }

    zero_array(prop->avg_data);
    zero_array(prop->std_dev_data);
    zero_array(prop->filter_fraction);
    clear_histogram(&prop->full_histogram);
    clear_histogram(&prop->filt_histogram);

    for (const auto p : ctx.properties) {
        if (p != prop && compare(p->name_buf, prop->name_buf)) {
            set_error_message(prop, "A property is already defined with that name!");
            prop->valid = false;
            return false;
        }
        if (p == prop) break;
    }

    if (!balanced_parentheses(prop->args_buf)) {
        set_error_message(prop, "Unbalanced parantheses!");
        prop->valid = false;
        return false;
    }

    DynamicArray<CString> args;

    // Extract big argument chunks
    const char* beg = prop->args_buf.beg();
    const char* end = prop->args_buf.beg();
    int count = 0;

    // Use space separation unless we are inside a parenthesis
    while (end != prop->args_buf.end()) {
        if (*end == '(')
            count++;
        else if (*end == ')')
            count--;
        else if (count == 0 && isspace(*end)) {
            CString arg = trim({beg, end});
            if (arg.size() > 0) {
                args.push_back(trim({beg, end}));
            }
            beg = end + 1;
        }
        end++;
    }
    if (beg != end) {
        CString arg = trim({beg, end});
        if (arg.size() > 0) {
            args.push_back(trim({beg, end}));
        }
    }

    if (args.size() == 0) {
        prop->valid = false;
        return false;
    }

    CString cmd = args[0];
    args = args.subarray(1);

    auto func = find_property_compute_func(cmd);
    if (!func) {
        set_error_message(prop, "Could not recognize command '%.*s'", cmd.length(), cmd.cstr());
        prop->valid = false;
        return false;
    }

    Range pre_range = prop->total_data_range;
    if (!func(prop, args, dynamic)) {
        prop->valid = false;
        return false;
    }

    if (prop->instance_data.count > 1) {
        if (prop->std_dev_data[0] == 0.f) {
            // Compute std dev of instances
            const float scl = 1.f / (float)prop->instance_data.count;
            for (int32 i = 0; i < num_frames; i++) {
                float sum = 0.f;
                for (const auto& inst : prop->instance_data) {
                    float x = inst.data[i] - prop->avg_data[i];
                    sum += x * x;
                }
                prop->std_dev_data[i] = math::sqrt(sum * scl);
            }
        }
    }

    if (pre_range != prop->total_data_range) {
        prop->filter = prop->total_data_range;
    }

    prop->full_histogram.value_range = prop->total_data_range;
    prop->filt_histogram.value_range = prop->total_data_range;
    prop->filter_dirty = true;
    prop->valid = true;

    ctx.current_property = nullptr;

    return true;
}

bool properties_dirty() {
    for (const auto p : ctx.properties) {
        if (p->data_dirty) return true;
        if (p->filter_dirty) return true;
    }

    return false;
}

void async_update(const MoleculeDynamic& dynamic, Range<int32> frame_filter, void (*on_finished)(void*), void* usr_data) {
    if (!dynamic) return;
    if (dynamic.trajectory.num_frames == 0) return;
    static Range<int32> prev_frame_filter{-1, -1};

    const bool frame_filter_changed = prev_frame_filter != frame_filter;
    const bool dirty_props = properties_dirty();

    if ((frame_filter_changed || dirty_props) && !ctx.thread_running && !ctx.stop_signal) {
        prev_frame_filter = frame_filter;
        ctx.thread_running = true;
        std::thread([&dynamic, frame_filter, on_finished, usr_data]() {
            Histogram tmp_hist;
            init_histogram(&tmp_hist, NUM_BINS);
            defer { free_histogram(&tmp_hist); };
            ctx.fraction_done = 0.f;

            // @NOTE IMPORTANT: This is the one 'true' frame count which should be used for properties.
            // It is important that this is used so that all properties data lengths are in sync
            // When dealing with dependencies.
            const int32 num_frames = dynamic.trajectory.num_frames;

            for (int32 i = 0; i < ctx.properties.size(); i++) {
                auto p = ctx.properties[i];
                ctx.fraction_done = (i / (float)ctx.properties.size());
                auto filter = p->filter;

                if (p->data_dirty) {
                    compute_property_data(p, dynamic, num_frames);

                    // recompute full histogram
                    clear_histogram(&p->full_histogram);
                    if (p->instance_data) {
                        for (const auto& inst : p->instance_data) {
                            compute_histogram(&p->full_histogram, inst.data);
                        }
                    } else {
                        compute_histogram(&p->full_histogram, p->avg_data);
                    }
                    normalize_histogram(&p->full_histogram, p->full_histogram.num_samples);

                    p->data_dirty = false;
                }

                if (ctx.stop_signal) break;

                if (p->filter_dirty) {
                    clear_histogram(&tmp_hist);
                    tmp_hist.value_range = p->filt_histogram.value_range;

                    int32 beg_idx = math::clamp((int32)frame_filter.x, 0, (int32)p->avg_data.size());
                    int32 end_idx = math::clamp((int32)frame_filter.y, 0, (int32)p->avg_data.size());

                    if (beg_idx != end_idx) {
                        // Since the data is probably showing, perform the operations on tmp data then copy the results
                        if (p->instance_data) {
                            for (const auto& inst : p->instance_data) {
                                compute_histogram(&tmp_hist, inst.data.subarray(beg_idx, end_idx - beg_idx));
                            }
                        } else {
                            compute_histogram(&tmp_hist, p->avg_data.subarray(beg_idx, end_idx - beg_idx));
                        }
                        normalize_histogram(&tmp_hist, tmp_hist.num_samples);
                        p->filt_histogram.bin_range = tmp_hist.bin_range;
                        p->filt_histogram.num_samples = tmp_hist.num_samples;
                        memcpy(p->filt_histogram.bins.ptr, tmp_hist.bins.ptr, p->filt_histogram.bins.size_in_bytes());
                    }

                    // Compute filter fractions for frames
                    if (p->instance_data) {
                        for (int32 j = 0; j < p->filter_fraction.size(); j++) {
                            float val = 0.f;
                            for (const auto& inst : p->instance_data) {
                                if (filter.x <= inst.data[j] && inst.data[j] <= filter.y) {
                                    val += 1.f;
                                }
                            }
                            p->filter_fraction[j] = val / (float)p->instance_data.count;
                        }
                    }

                    if (p->filter == filter) p->filter_dirty = false;
                }

                if (ctx.stop_signal) break;
            }

            if (on_finished) {
                on_finished(usr_data);
            }

            ctx.fraction_done = 1.f;
            ctx.thread_running = false;
            ctx.stop_signal = false;
        }).detach();
    }
}

bool thread_running() { return ctx.thread_running; }

void signal_stop() { ctx.stop_signal = true; }

void signal_stop_and_wait() {
    ctx.stop_signal = true;
    while (ctx.thread_running) {
        // possibly sleep 1 ms or so
    }
    ctx.stop_signal = false;
}

void signal_start() { ctx.stop_signal = false; }

float fraction_done() { return ctx.fraction_done; }

void visualize(const MoleculeDynamic& dynamic) {
    if (thread_running()) return;
    for (auto p : ctx.properties) {
        if (!p->enable_visualization) continue;
        CString cmd = extract_command(p->args_buf);
        auto entry = find_property_func_entry(COMPUTE_ID(cmd));
        if (entry && entry->visualize_func) {
            entry->visualize_func(*p, dynamic);
        }
    }
}

const Volume& get_density_volume() { return ctx.volume; }

VisualizationStyle* get_style() { return &ctx.style; }

Property* create_property(CString name, CString args) {
    Property* prop = (Property*)MALLOC(sizeof(Property));
    new (prop) Property();
    prop->name_buf = name;
    prop->args_buf = args;
    prop->valid = false;
    prop->data_dirty = true;
    prop->filter_dirty = true;

    init_histogram(&prop->full_histogram, NUM_BINS);
    clear_histogram(&prop->full_histogram);

    init_histogram(&prop->filt_histogram, NUM_BINS);
    clear_histogram(&prop->filt_histogram);

    ctx.properties.push_back(prop);
    return prop;
}

static void free_property_data(Property* prop) {
    signal_stop_and_wait();
    ASSERT(prop);
    free_histogram(&prop->full_histogram);
    free_histogram(&prop->filt_histogram);
    free_structure_data(&prop->structure_data);
    free_instance_data(&prop->instance_data);
    ctx.stop_signal = false;
}

void remove_property(Property* prop) {
    ASSERT(prop);
    signal_stop_and_wait();
    for (auto& p : ctx.properties) {
        if (p == prop) {
            free_property_data(prop);
            ctx.properties.remove(&p);
        }
    }
    ctx.stop_signal = false;
}

void remove_all_properties() {
    signal_stop_and_wait();
    for (auto prop : ctx.properties) {
        free_property_data(prop);
    }
    ctx.properties.clear();
    ctx.stop_signal = false;
}

void move_property_up(Property* prop) {
    if (ctx.properties.size() <= 1) return;
    signal_stop_and_wait();
    for (int32 i = 1; i < (int32)ctx.properties.size(); i++) {
        if (ctx.properties[i] == prop) {
            // swap
            Property* tmp = ctx.properties[i - 1];
            ctx.properties[i - 1] = ctx.properties[i];
            ctx.properties[i] = tmp;
            break;
        }
    }
    ctx.stop_signal = false;
}

void move_property_down(Property* prop) {
    if (ctx.properties.size() <= 1) return;
    signal_stop_and_wait();
    for (int32 i = 0; i < (int32)ctx.properties.size() - 1; i++) {
        if (ctx.properties[i] == prop) {
            // swap
            Property* tmp = ctx.properties[i + 1];
            ctx.properties[i + 1] = ctx.properties[i];
            ctx.properties[i] = tmp;
            break;
        }
    }
    ctx.stop_signal = false;
}

Array<Property*> get_properties() { return ctx.properties; }

Property* find_property(CString name) {
    for (auto p : ctx.properties) {
        if (compare(p->name_buf, name)) return p;
    }
    return nullptr;
}

void clear_property(Property* prop) {
    signal_stop_and_wait();
    ASSERT(prop);
    for (auto p : ctx.properties) {
        if (p == prop) {
            free_structure_data(&prop->structure_data);
            free_instance_data(&prop->instance_data);
            clear_histogram(&prop->full_histogram);
            prop->avg_data.clear();
        }
    }
    ctx.stop_signal = false;
}

void clear_all_properties() {
    signal_stop_and_wait();
    for (auto p : ctx.properties) {
        free_structure_data(&p->structure_data);
        free_instance_data(&p->instance_data);
        clear_histogram(&p->full_histogram);
        p->avg_data.clear();
    }
    ctx.stop_signal = false;
}

void set_all_property_flags(bool data_dirty, bool filter_dirty) {
    for (auto p : ctx.properties) {
        p->data_dirty |= data_dirty;
        p->filter_dirty |= filter_dirty;
    }
}

}  // namespace prop
