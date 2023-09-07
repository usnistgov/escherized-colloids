# Elements from Table 5.2.3 in "Tilings and Patterns" by Gruenbaum and Shephard.

import pandas as pd
import numpy as np

def get_factors(x):
    f = []
    for i in range(1, x + 1):
        if x % i == 0:
            f.append(i)

    return f[:-1] # Exclude trivial factors of 1 and x

def get_subgroups(group):
    subgroups = []
    if (group == 'c1'): # Identity
        return []
    elif (group.startswith("d")): # Dihedral group of order 2n (dn) 
        if group.endswith("inf"): # dinf = group of isotropic circle
            # Contains c_inf and other c, also d
            # However, for the purposes of this code, we can just return all the (minimal)
            # groups found in the table, which should suffice for all logic needed.
            return ['c1', 'c2', 'c3', 'c4', 'c6', 'd1', 'd2', 'd3', 'd4', 'd6']
        else:
            n = int(group.split("d")[1])
            subgroups += ["c{}".format(i) for i in get_factors(n)] + ["c{}".format(n)]
            subgroups += ["d{}".format(i) for i in get_factors(n)]
    elif (group.startswith("c")): # Cyclic group of order n (cn)
        n = int(group.split("c")[1]) 
        subgroups += ["c{}".format(i) for i in get_factors(n)]

    return subgroups

row_1 = [
    'PP1',
    'p1',
    'c1',
    ['Primitive'],
    'c2',
    [1, 41]
]

row_2 = [
    'PP2',
    'pg',
    'c1',
    ['Primitive'],
    'dinf',
    [2, 3, 43, 44]
]

row_3 = [
    'PP3',
    'pm',
    'c1',
    ['Primitive'],
    'dinf',
    [42]
]

row_4 = [
    'PP4',
    'pm',
    'd1',
    ['p1', 'pg', 'cm', '*'],
    'd2',
    [64]
]

row_5 = [
    'PP5',
    'cm',
    'c1',
    ['Primitive'],
    'dinf',
    [22, 45, 83]
]

row_6 = [
    'PP6',
    'cm',
    'd1',
    ['p1', 'pg'],
    'd2',
    [12, 14, 68]
]

row_7 = [
    'PP7',
    'p2',
    'c1',
    ['Primitive'],
    '',
    [4, 23, 46, 47, 84]
]

row_8 = [
    'PP8',
    'p2',
    'c2', 
    ['p1', '*'],
    '',
    [8, 57]
]

row_9 = [
    'PP9',
    'pgg',
    'c1',
    ['Primitive'],
    '',
    [5, 6, 25, 27, 51, 52, 53, 86]
]

row_10 = [
    'PP10',
    'pgg',
    'c2',
    ['pg'],
    'dinf',
    [9, 59]
]

row_11 = [
    'PP11',
    'pmg',
    'c1',
    ['Primitive'],
    '',
    [24, 49, 50, 85]
]

row_12 = [
    'PP12',
    'pmg', 
    'c2',
    ['pg', 'pm', 'pgg', '*'],
    'dinf',
    [58]
]

row_13 = [
    'PP13',
    'pmg', 
    'd1',
    ['pg', 'p2', 'pgg'],
    '',
    [13, 15, 66, 69]
]

row_14 = [
    'PP14',
    'pmm',
    'c1', 
    ['Primitive'],
    '',
    [48]
]

row_15 = [
    'PP15',
    'pmm',
    'd1', 
    ['pm', 'p2', 'pmg(2)', 'cmm', '*'],
    '',
    [65]
]

row_16 = [
    'PP16',
    'pmm',
    'd2', 
    ['p1', 'pg', 'pm(2)', 'cm', 'p2(3)', 'pgg', 'pmg(2)', 'cmm(3)', '*(2)'],
    '',
    [72]
]

row_17 = [
    'PP17',
    'cmm',
    'c1', 
    ['Primitive'],
    '',
    [54, 78]
]

row_18 = [
    'PP18',
    'cmm',
    'c2', 
    ['cm', 'pgg', 'pmm'],
    'dinf',
    [60]
]

row_19 = [
    'PP19',
    'cmm',
    'd1', 
    ['cm', 'p2', 'pgg', 'pmg'],
    '',
    [26, 67, 91]
]

row_20 = [
    'PP20',
    'cmm',
    'd2', 
    ['p1', 'pg', 'cm', 'p2(2)', 'pgg(2)', 'pmg'],
    '',
    [17, 74]
]

row_21 = [
    'PP21',
    'p3',
    'c1', 
    ['Primitive'],
    '',
    [7, 33]
]

row_22 = [
    'PP22',
    'p3',
    'c3', 
    ['p1', '*'],
    'c6',
    [10]
]

row_23 = [
    'PP23',
    'p31m',
    'c1', 
    ['Primitive'],
    '',
    [30, 38]
]

row_24 = [
    'PP24',
    'p31m',
    'c3', 
    ['cm', 'p3m1'],
    'dinf',
    [89]
]

row_25 = [
    'PP25',
    'p31m',
    'd1', 
    ['p3'],
    '',
    [16, 36]
]

row_26 = [
    'PP26',
    'p31m',
    'd3', 
    ['p1', 'pg', 'cm', 'p3(2)'],
    'd6',
    [18]
]

row_27 = [
    'PP27',
    'p3m1',
    'c1', 
    ['Primitive'],
    '',
    [87]
]

row_28 = [
    'PP28',
    'p3m1',
    'd1', 
    ['p3'],
    '',
    [35]
]

row_29 = [
    'PP29',
    'p3m1',
    'd3', 
    ['p1', 'pg', 'cm', 'p3(2)', 'p31m'],
    'd6',
    [19]
]

row_30 = [
    'PP30',
    'p4',
    'c1', 
    ['Primitive'],
    '',
    [28, 55, 79]
]

row_31 = [
    'PP31',
    'p4',
    'c2', 
    ['*'],
    'c4',
    [61]
]

row_32 = [
    'PP32',
    'p4',
    'c4', 
    ['p1', 'p2(3)', '*(2)'],
    'dinf',
    [62]
]

row_33 = [
    'PP33',
    'p4g',
    'c1', 
    ['Primitive'],
    '',
    [56, 81]
]

row_34 = [
    'PP34',
    'p4g',
    'c4', 
    ['pg', 'cm', 'pgg(2)', 'pmm', 'cmm'],
    'dinf',
    [63]
]

row_35 = [
    'PP35',
    'p4g',
    'd1', 
    ['pgg', 'p4'],
    '',
    [29, 71]
]

row_36 = [
    'PP36',
    'p4g',
    'd2', 
    ['pg', 'pgg', 'p4(2)'],
    'd4',
    [73]
]

row_37 = [
    'PP37',
    'p4m',
    'c1', 
    ['Primitive'],
    '',
    [80]
]

row_38 = [
    'PP38',
    'p4m',
    'd1', 
    ['cmm', 'p4', 'p4g', '*'],
    '',
    [82]
]

row_39 = [
    'PP39',
    'p4m',
    'd1', 
    ['pmm', 'p4', 'p4g'],
    '',
    [70]
]

row_40 = [
    'PP40',
    'p4m',
    'd2', 
    ['cm', 'pgg', 'pmm', 'cmm', 'p4(2)', 'p4g(2)', '*'],
    'd4',
    [75]
]

row_41 = [
    'PP41',
    'p4m',
    'd4', 
    ['p1', 'pg(2)', 'pm(2)', 'cm(2)', 'p2(3)', 'pgg(3)', 
     'pmg(3)', 'pmm(3)', 'cmm(4)', 'p4(3)', 'p4g(3)', '*(2)'],
    '',
    [76]
]

row_42 = [
    'PP42',
    'p6',
    'c1', 
    ['Primitive'],
    '',
    [21, 31, 39, 88]
]

row_43 = [
    'PP43',
    'p6',
    'c2', 
    ['p3'],
    'dinf',
    [34]
]

row_44 = [
    'PP44',
    'p6',
    'c3', 
    ['p2', '*'],
    'dinf',
    [90]
]

row_45 = [
    'PP45',
    'p6',
    'c6', 
    ['p1', 'p2(2)', 'p3(2)'],
    'dinf',
    [11]
]

row_46 = [
    'PP46',
    'p6m',
    'c1', 
    ['Primitive'],
    '',
    [77]
]

row_47 = [
    'PP47',
    'p6m',
    'd1', 
    ['p3m1', 'p6'],
    '',
    [92]
]

row_48 = [
    'PP48',
    'p6m',
    'd1', 
    ['p31m', 'p6'],
    '',
    [32, 40]
]

row_49 = [
    'PP49',
    'p6m',
    'd2', 
    ['p3', 'p31m', 'p3m1', 'p6'],
    '',
    [37]
]

row_50 = [
    'PP50',
    'p6m',
    'd3', 
    ['cm', 'pgg', 'pmg', 'cmm', 'p2', 'p31m', 'p3m1', 'p6(2)'],
    '',
    [93]
]

row_51 = [
    'PP51',
    'p6m',
    'd6', 
    ['p1', 'pg(2)', 'cm(2)', 'p2(2)', 'pgg(3)', 'pmg(2)', 'cmm', 'p3(2)', 'p31m(2)', 'p3m1', 'p6'],
    '',
    [20]
]

col_names = ['Pattern Type', 
             'Symmetry Group '+r'$S(\mathcal{P})$', 
             'Induced Group, '+r'$S(\mathcal{P}|M)$',
             'Motif-transitive Proper Subgroups of '+r'$S(\mathcal{P})$',
             'Minimal Forbidden Supergroups',
             'Isohedral Tiling Type IH'
            ]

symmetry_table = pd.DataFrame(
    data=[row_1, row_2, row_3, row_4, row_5, row_6, row_7, row_8, row_9, row_10,
          row_11, row_12, row_13, row_14, row_15, row_16, row_17, row_18, row_19, row_20,
          row_21, row_22, row_23, row_24, row_25, row_26, row_27, row_28, row_29, row_30,
          row_31, row_32, row_33, row_34, row_35, row_36, row_37, row_38, row_39, row_40,
          row_41, row_42, row_43, row_44, row_45, row_46, row_47, row_48, row_49, row_50,
          row_51
         ],
    columns=col_names
)

def prioritize(motif_point_symmetry):
    mask1 = np.array(symmetry_table[col_names[2]] == motif_point_symmetry)
    safe = symmetry_table[mask1]
    
    # Induced group is subgroup of S(M) - if it is equal then it is safe instead
    mask2 = np.array([x in get_subgroups(motif_point_symmetry) for x in symmetry_table[col_names[2]]])
    # S(M) cannot be forbidden supergroup nor have a subgroup that is
    mask3 = np.array([x in [motif_point_symmetry] + get_subgroups(motif_point_symmetry) for x in symmetry_table[col_names[4]]])
    dangerous = symmetry_table[mask2 & ~mask3]
    
    # Either the induced group is not a subgroup of the motif, S(M) - if equal it is safe
    # or if it is, then the minimal forbidden supergroup belongs to the motif group.
    # This is the opposite of the dangerous mask and exclude safe ones, so they are the remainder
    strictly_forbidden = symmetry_table[~mask1 & ~(mask2 & ~mask3)]
    
    return safe, dangerous, strictly_forbidden

def watch_out_for(ih_type, dangerous, motif_point_symmetry):
    """
    For "dangerous" cases, the motif has more symmetry than we are intentionally inducing.
    The code will place the motif so it induces at least what is required by the tile's
    properties (e.g., on a mirror line).  "Forbidden" supergroups are just groups that are
    always strictly forbidden, but "dangerous" ones may work, though if the right
    set of conditions (tile shape, etc.) are met you might actually end up having (some) of 
    those extra symmetries in the motif become part of the pattern (or other effects which
    result in additional global symmetries). In fact, this might only change the henomeric 
    type without changing the symmetry group (see G&S) but we will disallow this
    to be more conservative.
    
    Logic:
    1. Start with the S(P1) you are targeting and would usually get.
    2. That S(P1) should be treated as the MTPS of another pattern, S(P2), in which the 
        induced group S(P2|M) is a proper subgroup (inclusive) of S(M). 
    - This is because you cannot allow S(P2) to require induction of more symmetry than the 
        motif has, since this is impossible for the motif provide. 
    3. We also require that S(P2|M) > S(P1|M), since we are considering the case where 
        the new pattern is exploiting more of the motif's "extra" symmetries than before.
        Note, S(M) > S(P1|M) - otherwise this would not be a dangerous case.
    4. We should also remove any cases where the minimal forbidden supergroup is a proper subgroup 
        of S(M) also since this would always forbid this motif.
        
    This will return rows in the table that tell you what symmetries you might end up with
    instead of what you were hoping for.
    """
    df = dangerous.loc[[ih_type in x for x in dangerous[col_names[5]].values]]
    if df.empty:
        raise ValueError("IH tile "+str(ih_type)+" is not a dangerous type")

    target_group = df[col_names[1]].values[0] # Wallpaper group originally targeted
    point_symm = df[col_names[2]].values[0] # The motif has more than this, otherwise it would have been considered "safe" by our algorithm
    
    def rename(values, group):
        # 1. Remove any parenthases
        new_v = []
        for v in values:
            new_v.append(v.split('(')[0])
            
        # 2. Explicitly rename isomorphic ones
        for i in range(len(new_v)):
            if new_v[i] == '*':
                new_v[i] = group
                
        return new_v

    possibilities = symmetry_table.loc[[target_group in rename(v,g) for g,v in symmetry_table[[col_names[1], col_names[3]]].values]]
    
    if (len(possibilities) > 0):
        # Remove any possibilities that induce more/different than motif has to begin with (2)
        mask1 = np.array([x in [motif_point_symmetry] + get_subgroups(motif_point_symmetry) 
                          for x in possibilities[col_names[2]]])
        # Must induce MORE than the point_symm, i.e. point_symm is a SUBgroup (3)
        mask2 = np.array([point_symm in get_subgroups(x) 
                          for x in possibilities[col_names[2]]])
        # Remove any forbidden cases (4)
        mask3 = np.array([x in [motif_point_symmetry] + get_subgroups(motif_point_symmetry) 
                          for x in possibilities[col_names[4]]])
        return possibilities[mask1 & mask2 & ~mask3]
    else:
        return possibilities # Just an empty table

if __name__ == "__main__":

    def sanity_checks(df):
        # 51 rows
        assert(df.shape[0] == 51)
    
        # Check all IH tiles used once and only once
        used = {}
        for x in df['Isohedral Tiling Type IH'].values:
            for v in x:
                assert(v not in used), '{}'.format(v)+ ' is repeated'
                used[v] = True
        assert(np.all(np.array(sorted(used.keys())) == np.arange(1, 93+1)))
    
        # Forbidden groups not equal induced
        for a, b in df[[col_names[2], col_names[4]]].values:
            assert(a != b)
        
        # All forbidden groups are supergroups of induced
        for a, b in df[[col_names[2], col_names[4]]].values:
            if b != '':
                assert(a in get_subgroups(b))
        
        # No MTPS are repeated
        for x in df[col_names[3]].values:
            _, c = np.unique(x, return_counts=True)
            assert(np.all(c == 1))

    sanity_checks(symmetry_table)

