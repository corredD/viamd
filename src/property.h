#pragma once

#include <core/common.h>
#include <core/types.h>
#include <core/array_types.h>
#include <core/string_types.h>
#include <core/bitfield.h>

struct MoleculeDynamic;


/*

Computation Commands:

angle(one, one, one)

dihedral(one, one, one, one)

distance(one, one)
min_distance(many, many)
max_distance(many, many)
mean_distance(many, many)

rmsd(many)

position(one)
mean_position(many)


Structure Matching Commands:
1,2,3,4 will match atomic indices within whatever context that is given
* will match all
1:4 will match 1 to 4
2:* will match 2 to N

residue(RES) will match all atoms which belong to any residue named RES
residue(1:4) will match atoms within residues 1 to 4
residue(ALA LYS GLY) will match all atoms which belong to any residues named ALA, LYS or GLY

chain(A) will match all atoms which belong to any chain named A
chain(2) will match all atoms which belong to chain with index 2


Examples:

a1 = angle(1, 2, 3)                                   // Angle between atoms 1 2 3

d1 = for_each_residue(RES, distance(com(res), res(1)))    // for each residue that matches calculate the distance between the com of the residue and its first atom
    d1[0]...
    d1[1]...
    d1[2]...
  
d2 = for_each_residue(RES ALN VYS, min_distance(atom(res, *), com(atom(res, *))) // OK: for each residue that matches RES it locally computes the minimal distance between all atoms and its center of mass
    d2[0]...
    d2[1]...
    d2[2]...

d3 = distance(residue(RES), 1)                         // ERROR! argument 1: 'residue(res)' matches many atoms
d4 = distance(com(residue(ALA LYS GLY)), 1)            // OK: computes distance between com of all residues matching ALA or LYS or GLY and atom 1

r1 = for_each_residue(1:4, rmsd(*))                  // computes rmsd for residues 1, 2, 3, 4
    r1[0]...
    r1[1]...
    r1[2]...

p1 = position(2)                                       // computes position of atom 2
    p1.x...
    p1.y...
    p1.z...

p2 = for_each_residue(ALA LYS GLY, position(com(res(*)))            // computes com of each residue and stores as x y z child properties
    p2.x[]...
        p2.x[0]...
        p2.x[1]...
        p2.x[2]...
    p2.y[]...
        p2.y[0]...
        p2.y[1]...
        p2.y[2]...
    p2.z[]...
        p2.z[0]...
        p2.z[1]...
        p2.z[2]...
*/

namespace prop {

/*
// HISTOGRAM
struct Histogram {
    Array<float> bins = {};
    Range<float> val_range = {0, 0};
    Range<float> bin_range = {0, 0};
    int32 num_samples = 0;
};

void init_histogram(Histogram* hist, int32 num_bins);
void free_histogram(Histogram* hist);

void compute_histogram(Histogram* hist, Array<const float> data);
void compute_histogram(Histogram* hist, Array<const float> data, Range<float> filter);

void clear_histogram(Histogram* hist);
void normalize_histogram(Histogram* hist);
*/

using ID = uint64;
constexpr ID INVALID_ID = 0;

ID        create_property(CString name, CString args, CString unit = "", bool is_periodic = false);
bool      remove_property(ID prop_id);
Array<ID> get_properties();

CString   get_name(ID prop_id);
CString   get_unit(ID prop_id);
bool      get_periodic(ID prop_id);
bool      get_valid(ID prop_id);
CString   get_error_msg(ID prop_id);
Bitfield  get_atoms(ID prop_id);
ID        get_parent(ID prop_id);
Array<const float> get_data(ID prop_id);
Array<const ID>    get_children(ID prop_id);

// Async functionality
void async_update(const MoleculeDynamic& dynamic, Range<int32> frame_filter = {0, 0}, void (*on_finished_callback)(void*) = nullptr, void* usr_data = nullptr);
bool thread_running();
void signal_stop();
void signal_stop_and_wait();
float fraction_done();

// Commands
Array<CString> get_property_commands();
Array<CString> get_structure_commands();

// Property Manager
void initialize();
void shutdown();

}  // namespace prop
