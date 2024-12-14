from sage.all import euler_phi, lazy_attribute, point, line, polygon, frac, floor, lcm, mod, cartesian_product_iterator, ZZ, QQ, PolynomialRing, OrderedPartitions, srange, prime_range, prime_pi, next_prime, previous_prime, gcd, conway_polynomial, GF, binomial, cached_function, infinity, latex
from sage.databases.cremona import cremona_letter_code
from lmfdb import db
from lmfdb.galois_groups.transitive_group import knowl_cache, transitive_group_display_knowl
from flask import url_for

from collections import defaultdict, Counter
import itertools
import re
import os
FAMILY_RE = re.compile(r'(\d+(?:\.\d+\.\d+\.\d+)?)-(?:(\d+)\.(\d+(?:_\d+)*))?')

columns = [
    ("label", "text"),
    ("ctr", "smallint"),
    ("base", "text"),
    ("base_aut", "smallint"),
    ("e0", "smallint"),
    ("e", "smallint"),
    ("e_absolute", "smallint"),
    ("f0", "smallint"),
    ("f", "smallint"),
    ("f_absolute", "smallint"),
    ("n0", "smallint"),
    ("n", "smallint"),
    ("n_absolute", "smallint"),
    ("w", "smallint"),
    ("p", "integer"),
    ("c0", "smallint"),
    ("c", "smallint"),
    ("c_absolute", "smallint"),
    ("rf0", "smallint[]"), # TODO
    ("visible", "text"),
    ("slopes", "text"),
    ("heights", "text"),
    ("scaled_heights", "text"),
    ("rams", "text"),
    ("scaled_rams", "text"),
    ("field_count", "integer"),
    ("packet_count", "integer"),
    ("ambiguity", "smallint"),
    ("mass", "double precision"),
    ("mass_stored", "text"),
    ("mass_missing", "double precision"),
    ("mass_display", "text"),
    ("all_stored", "boolean"),
    ("slope_multiplicities", "smallint[]"),
    ("wild_segments", "smallint"),
    ("poly", "text"),
    ("label_absolute", "text"),
    ("rams_absolute", "text"),
    ("top_slope", "double precision"),
]

cache_cols = ["p", "e", "f", "c", "family", "packet", "label", "new_label", "coeffs", "galT", "galois_label", "galois_degree", "slopes", "ind_of_insep", "associated_inertia", "t", "u", "aut", "visible", "hidden", "subfield", "residual_polynomials", "rf"]
def cache_key(rec):
    return rec["p"], rec["f"], rec["e"], rec["visible"]

def get_new_labels():
    # The order within a subfamily is currently wrong; we will need to rerun this once we have polredpadics
    R = PolynomialRing(ZZ, ["t", "z"])
    by_fam = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
    by_pack = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(list)))))
    for rec in db.lf_fields.search(
            {},
            ["label", "p", "f", "e", "c", "visible", "residual_polynomials", "galT", "t", "u", "hidden"]):
        p, f, e, c, galT = rec["p"], rec["f"], rec["e"], rec["c"], rec["galT"]
        if rec["visible"] == "[]":
            vis = ()
        else:
            vis = tuple(QQ(x) for x in rec["visible"][1:-1].split(", "))
        if rec["hidden"] == "[]":
            hid = ()
        else:
            hid = tuple(QQ(x) for x in rec["hidden"][1:-1].split(", "))
        rpol = tuple(R(pol) for pol in rec["residual_polynomials"])
        by_fam[p,f,e,c][vis][rpol].append(rec["label"])
        by_pack[p,f,e,c][vis][rpol][galT][hid].append(rec["label"])
    new_label = {}
    for (p,f,e,c), A in by_fam.items():
        fam_ctr = 0
        for vis in sorted(A):
            fam_code = cremona_letter_code(fam_ctr)
            fam_ctr += 1
            B = A[vis]
            sub_ctr = 0
            for rpol in sorted(B):
                sub_ctr += 1
                pack_lookup = {}
                for galT, C in by_pack[p,f,e,c][vis][rpol].items():
                    # NOTE: hctr could change if we are missing fields in this subfamily
                    for hctr, hid in enumerate(sorted(C)):
                        for label in C[hid]:
                            pack_lookup[label] = f"{galT}{cremona_letter_code(hctr)}"
                for i, label in enumerate(B[rpol], 1):
                    family = f"{p}.{f}.{e}.{c}{fam_code}"
                    nl = family + f"{sub_ctr}.{i}"
                    packet = family + f"{sub_ctr}.{pack_lookup[label]}"
                    new_label[label] = (nl, family, packet, fam_ctr, sub_ctr, i)
    return new_label

def write_new_labels(fname):
    with open(fname, "w") as F:
        _ = F.write("label|new_label|family|packet|ctr_family|ctr_subfamily|ctr\ntext|text|text|text|integer|integer|integer\n\n")
        for label, (new_label, family, packet, fam_ctr, sub_ctr, ctr) in get_new_labels().items():
            _ = F.write(f"{label}|{new_label}|{family}|{packet}|{fam_ctr}|{sub_ctr}|{ctr}\n")

def str_to_QQtup(s):
    if s == "[]":
        return ()
    return tuple(QQ(x) for x in s[1:-1].split(", "))

def get_filtration(rec, vis, cans):
    label =rec["new_label"]
    p, _ = label.split(".", 1)
    p = ZZ(p)
    f, e = ZZ(rec["f"]), ZZ(rec["e"])
    w, etame = e.val_unit(p)
    filt = [f"{p}.{f}.1.0a1.1"]
    poss = defaultdict(list)
    for K in rec["subfield"]:
        _, f0, e0, _, _ = K.split(".")
        f0, e0 = ZZ(f0), ZZ(e0)
        if f0 == f and e0 % etame == 0:
            poss[f, e0].append(K)
    if etame == 1:
        filt.append(filt[0])
    elif w == 0:
        filt.append(label)
    else:
        assert len(poss[f, etame]) == 1
        filt.append(poss[f, etame][0])
    ctr = 1
    for i, s in enumerate(vis[:-1], 1):
        if s == vis[i]:
            ctr += 1
        else:
            e1 = etame * p**i
            for K in poss[f, e1]:
                if K in cans[(p, f, e1) + vis[:i]]:
                    filt.extend([K] * ctr)
                    break
            else:
                raise RuntimeError("Canonical field not found")
            ctr = 1
    filt.extend([label] * (ctr - 1))
    return filt

def write_filtration(fname):
    old_to_new = {rec["label"]: rec["new_label"] for rec in db.lf_fields.search({}, ["label", "new_label"])}
    field_cache = defaultdict(list)
    for rec in db.lf_fields.search(
            {},
            ["family", "packet", "label", "new_label", "coeffs", "galT", "galois_label", "galois_degree", "slopes", "ind_of_insep", "associated_inertia", "t", "u", "aut", "hidden", "c"]):
        field_cache[rec["family"]].append(rec)
    cans = defaultdict(set)
    filts = {}
    for rec in db.lf_fields.search({}, ["p", "f", "e", "visible", "new_label", "subfield"]):
        label = rec["new_label"]
        p, f, e = ZZ(rec["p"]), ZZ(rec["f"]), ZZ(rec["e"])
        w, etame = e.val_unit(p)
        vis = str_to_QQtup(rec["visible"])
        rec["subfield"] = [old_to_new[sub] for sub in rec["subfield"]]
        cans[(p, f, e) + vis].add(label)
        filts[label] = get_filtration(rec, vis, cans)
    with open(fname, "w") as F:
        _ = F.write("new_label|canonical_filtration\ntext|text[]\n\n")
        for label, filt in filts.items():
            _ = F.write(f"{label}|{{{','.join(filt)}}}\n")


    #tame = []
    #for rec in db.lf_fields.search({"n": 1}, ["p", "new_label"]):
    #    old_to_new[str(rec["p"])] = rec["new_label"]
    #tame_sub = {a: old_to_new[b.split("-")[0]] for (a,b) in db._execute(SQL("SELECT label, family FROM lf_fields_old1"))}
    #tame_aut = {rec["new_label"]: rec["aut"] for rec in db.lf_fields.search({}, ["new_label", "aut"])}
    #wild_families = defaultdict(list)
    #tame_e = {}
    #tame_f = {}
    #seen = set()
    #for rec in db.lf_fields.search({}, ["label", "family", "p", "e", "f", "aut", "visible"]):
    #    family = rec["family"]
    #    p, f, e = ZZ(rec["p"]), ZZ(rec["f"]), ZZ(rec["e"])
    #    w, etame = e.val_unit(p)
    #    if family not in seen:
    #        if rec["visible"] == "[]":
    #            slopes = []
    #        else:
    #            slopes = [(QQ(x)-1)*etame for x in rec["visible"][1:-1].split(", ")]
    #        tame = tame_sub[rec["label"]]
    #        base_aut = tame_aut[tame]
    #        yield pAdicSlopeFamily(tame, e=etame, f=f, base_aut=base_aut, slopes=slopes, label=family, field_cache=field_cache)
    #        seen.add(family)

    #by_base = defaultdict(set)
    #for x in wild_families:
    #    base, slopes = x.split("-")
    #    by_base[base].add(slopes.count("_") + 1)
    #wild = []
    #for base, Ws in by_base.items():
    #    for w in sorted(Ws):
    #        yield from pAdicSlopeFamily.families(base, w, tame_e=tame_e[base], tame_f=tame_f[base], tame_aut=tame_aut[base], field_cache=field_cache)

def get_basepairs():
    raise RuntimeError("Old version")
    basepairs = set()
    def getn(label):
        _, e, f, _ = label.split(".", 3)
        return ZZ(e) * ZZ(f)
    for rec in db.lf_fields.search({}, ["new_label", "canonical_filtration", "p", "n"]):
        p, n = rec["p"], rec["n"]
        basepairs.add((p, n))
        for base in rec["canonical_filtration"]:
            m = getn(base)
            if m not in [1,n]:
                basepairs.add((base, n//m))
    return basepairs

def absolute_basepairs():
    basepairs = []
    for p in prime_range(200):
        for n in range(1,48):
            basepairs.append((p,n))
    return basepairs

def relative_basepairs():
    basepairs = []
    for rec in db.lf_fields.search({"n": {"$gt": 1, "$lt": 16}}, ["new_label", "n", "p", "e"]):
        p, n, e = ZZ(rec["p"]), ZZ(rec["n"]), ZZ(rec["e"])
        for m in srange(2, 48//n):
            # We used to only support extensions where the base was part of the canonical filtration of the absolute field.
            #if e == 1 or m.is_power_of(p):
            # If the base is not unramified, we can only support totally wildly ramified degrees
            basepairs.append((rec["new_label"], m))
    return basepairs

def get_field_data():
    field_data = {rec["new_label"]: rec for rec in db.lf_fields.search({}, cache_cols) if "new_label" in rec}
    for rec in field_data.values():
        rec["visible"] = str_to_QQtup(rec["visible"])
    return field_data

def get_caches():
    field_cache = defaultdict(list)
    field_data = {}
    for i, rec in enumerate(db.lf_fields.search({}, cache_cols)):
        if "new_label" not in rec:
            continue
        if i and i % 10000 == 0: print(i)
        field_cache[cache_key(rec)].append(rec)
        field_data[rec["new_label"]] = dict(rec)
        field_data[rec["new_label"]]["visible"] = str_to_QQtup(rec["visible"])
    return field_cache, field_data

def get_db_families(basepairs, field_cache, field_data):
    for i, (base, n) in enumerate(basepairs):
        if i and i % 100 == 0:
            print(i, base, n)
        yield from pAdicSlopeFamily.families(base, n, field_cache=field_cache, field_data=field_data)

def get_polynomials(basepairs, field_cache, field_data):
    for fam in get_db_families(basepairs, field_cache, field_data):
        yield from fam

def write_col(family, col, typ):
    out = getattr(family, col)
    if out is None:
        return r"\N"
    if typ == "double precision":
        out = float(out)
    out = str(out)
    if typ.endswith("[]"):
        out = "{" + out[1:-1] + "}"
    return out

def write_db_families(famfile, labelfile, corfile, basepairs, field_cache=None, field_data=None):
    # All absolute basepairs should come before relative basepairs, and the absolute basepair corresponding to every relative basepair present must be included.
    if field_cache is None or field_data is None:
        field_cache, field_data = get_caches()
    Ztz = PolynomialRing(PolynomialRing(ZZ, "t"), "z")
    absolute_lookup = {}
    with open(famfile, "w") as Ffam:
        _ = Ffam.write("|".join(col for (col,typ) in columns) + "\n")
        _ = Ffam.write("|".join(typ for (col,typ) in columns) + "\n\n")
        with open(labelfile, "w") as Flab:
            _ = Flab.write("label|family|subfamily|new_label\ntext|text|text|text\n\n")
            with open(corfile, "w") as Fcor:
                _ = Fcor.write("family|field\ntext|text\n\n")
                for family in get_db_families(basepairs, field_cache, field_data):
                    key = (family.p, family.f_absolute, family.e_absolute) + tuple(family.visible)
                    if family.n0 == 1:
                        #fields = sorted(family.fields, key=family.field_sort_key)
                        fields = family.fields
                        respoly = sorted((tuple(Ztz(f) for f in rec["residual_polynomials"]), i) for (i, rec) in enumerate(fields))
                        subctr = 0
                        sublook = {}
                        prev = None
                        for ftup, i in respoly:
                            rec = fields[i]
                            if ftup != prev:
                                prev = ftup
                                sublook[ftup] = subctr
                                subctr += 1
                                ctr = 1
                            rec["subfamily"] = f"{family.label}{subctr}"
                            rec["new_label"] = f"{family.label}{subctr}.{ctr}"
                            ctr += 1
                        for rec in fields:
                            _ = Flab.write(f"{rec['label']}|{family.label}|{rec['subfamily']}|{rec['new_label']}\n")
                        family.label_absolute = absolute_lookup[key] = family.label
                    else:
                        family.label_absolute = absolute_lookup[key]
                    for rec in family.fields:
                        _ = Fcor.write(f"{family.label}|{rec['new_label']}\n")
                    _ = Ffam.write("|".join(write_col(family, col, typ) for (col,typ) in columns) + "\n")

def write_polynomials(poly_folder, basepairs, field_cache, field_data, batch_size=64):
    for fam in get_db_families(basepairs, field_cache, field_data):
        ctr = 0
        F = None
        try:
            for poly in fam:
                if ctr % batch_size == 0:
                    if F is not None:
                        F.close()
                        F = None
                    F = open(poly_folder / f"{fam.label}.{ctr // batch_size}", "w")
                ctr += 1
                _ = F.write(f"{str(poly).replace(' ','')}\n")
        finally:
            if F is not None:
                F.close()

def rams_to_heights(rams, p):
    w = len(rams)
    return [(p-1)*sum(p**(k-j) * rams[j] for j in range(k+1)) for k in range(w)]

def slopes_to_heights(slopes, p, etame):
    heights = []
    h = 0
    phipk = etame * (p - 1)
    for s in slopes:
        h += phipk * s
        heights.append(h)
        phipk *= p
    return heights

def heights_to_rams(heights, p):
    w = len(heights)
    return [heights[0] / (p - 1)] + [(heights[k] - p*heights[k-1]) / (p - 1) for k in range(1,w)]

def heights_to_slopes(heights, p, etame):
    w = len(heights)
    return [heights[0] / (etame * (p-1))] + [(heights[k] - heights[k-1]) / (etame * euler_phi(p**(k+1))) for k in range(1,w)]

def sort_rams(rams, p):
    # Given rams that are not sorted, return the sorted version.  See David Robert's notes, section 2.8
    # We just use a simple bubble sort since the number of rams will always be very small (at most 5 for (2,32))
    rams = list(rams) # Shouldn't modify input
    w = len(rams)
    for i in range(1,w):
        for j in range(w-i):
            a, b = rams[j:j+2]
            if a > b:
                rams[j:j+2] = [b, b + p*(a-b)]
    return rams

def check_conductor_formula(rec, efcheight_cache):
    raise RuntimeError("Old version")
    def ef(label):
        _, e, f, _ = label.split(".", 3)
        return ZZ(e), ZZ(f)
    bases = [x for x in set(rec["canonical_filtration"]) if ef(x) != (1,1) and x != rec["new_label"]]
    p = rec["p"]
    e, f, c, height = efcheight_cache[rec["new_label"]]
    n = e * f
    h = height[-1] if height else 0
    w, etame = e.val_unit(p)
    assert c == f * (h + e - 1)
    for base in bases:
        e0, f0, c0, height0 = efcheight_cache[base]
        n0 = e0 * f0
        h0 = height0[-1] if height0 else 0
        w0, e0tame = e0.val_unit(p)
        assert c0 == f0 * (h0 + e0 - 1)
        if etame == e0tame:
            # relative extension is totally ramified
            assert c == c0 * (n // n0) + (h - h0)
        else:
            # unsure here
            assert c == c0 * (n // n0) + (h - h0)

def check_all_conductors():
    def get_heights(rec):
        slopes = [x - 1 for x in str_to_QQtup(rec["visible"])]
        p = ZZ(rec["p"])
        w, etame = ZZ(rec["e"]).val_unit(p)
        return slopes_to_heights(slopes, p, etame)
    efcheight_cache = {rec["new_label"]: (ZZ(rec["e"]), ZZ(rec["f"]), ZZ(rec["c"]), get_heights(rec)) for rec in db.lf_fields.search({}, ["new_label", "p", "e", "f", "c", "visible"])}
    for rec in db.lf_fields.search({}, ["new_label", "canonical_filtration", "p"]):
        check_conductor_formula(rec, efcheight_cache)

def check_respoly(basepairs, field_cache, field_data):
    for fam in get_db_families(basepairs, field_cache, field_data):
        if fam.n0 == 1 and fam.f == 1 and fam.etame == 1 and fam.n < 16:
            kz = GF(fam.p)['z']
            assert set(tuple(kz(f) for f in rec["residual_polynomials"]) for rec in fam.fields) == set(tuple(rp) for rp in fam.respoly_iter())

def check_respoly2(basepairs, field_cache, field_data):
    for fam in get_db_families(basepairs, field_cache, field_data):
        if fam.n0 == 1 and fam.f == 1 and fam.n < 16:
            kz = GF(fam.p)['z']
            assert set(tuple(kz(f) for f in rec["residual_polynomials"]) for rec in fam.fields) == set(tuple(rp) for rp in fam.respoly2_fixed_iter())


def respoly_class(rp, slopes):
    k = rp[0].base_ring()
    z = rp[0].parent().gen()
    cls = set()
    for delta in k:
        if delta == 0:
            continue
        gamma = 1
        rvec = []
        for s, rpoly in zip(slopes, rp):
            h, e = s.numerator(), s.denominator()
            gamma *= delta**(-h * rpoly.degree())
            rvec.append(gamma * rpoly(delta**h * z))
        cls.add(tuple(rvec))
    return cls

@cached_function
def respoly_rep(rp, slopes):
    print("Computing")
    cls = respoly_class(rp, slopes)
    rep = min(cls)
    for rp1 in cls:
        respoly_rep.set_cache(rep, (rp1, slopes))
    return rep

class pAdicSlopeFamily:
    def __init__(self, base, f=1, etame=1, slopes=[], heights=[], rams=[], base_data=None, field_cache=None):
        if base_data is None:
            base_data = db.lf_fields.lucky({"new_label": base}, cache_cols)
        self.base, self.f, self.etame, self.base_aut, self.oldbase = base, f, etame, base_data["aut"], base_data["label"]
        p, f0, e0, base_fam, base_i = base.split(".")
        p, f0, e0 = ZZ(p), ZZ(f0), ZZ(e0)
        w0, e0tame = e0.val_unit(p)
        self.c0 = c0 = ZZ(re.split(r"\D", base_fam)[0])
        self.rf0 = base_data["rf"]
        # For now, these slopes are Serre-Swan slopes, not Artin-Fontaine slopes
        assert p.is_prime()
        self.w = w = max(len(L) for L in [slopes, heights, rams])
        # We support tamely ramified fields by specifying empty slopes/rams/heights
        # slopes/rams -> heights -> rams/slopes
        if rams:
            heights = rams_to_heights(rams, p)
        if slopes:
            heights = slopes_to_heights(slopes, p, etame)
        if w and not rams:
            rams = heights_to_rams(heights, p)
        if w and not slopes:
            slopes = heights_to_slopes(heights, p, etame)
        self.slopes = slopes
        self.artin_slopes = [s + 1 for s in slopes]
        # Visible is absolute, for use in matching with fields
        base_rams = base_data.get("rams")
        if base_rams is None:
            if e0.gcd(p) == 1:
                base_rams = []
            else:
                base_slopes = [x - 1 for x in str_to_QQtup(db.lf_fields.lucky({"new_label":base}, "visible"))]
                base_rams = heights_to_rams(slopes_to_heights(base_slopes, p, 1), p)
        if base_rams:
            self.rams_absolute = sort_rams([etame * ram for ram in base_rams] + rams, p)
            self.visible = [x + 1 for x in heights_to_slopes(rams_to_heights(self.rams_absolute, p), p, e0tame*etame)]
        else:
            self.rams_absolute = rams
            self.visible = [x/e0tame + 1 for x in slopes]
        self.heights = heights
        self.rams = rams
        self.pw = pw = p**w
        self.e = e = etame * pw
        self.n = n = e * f
        self.e0 = e0
        self.e_absolute = e * e0
        self.f0 = f0
        self.f_absolute = f * f0
        self.n0 = n0 = e0 * f0
        self.n_absolute = n * n0
        self.p = p
        self.field_cache = field_cache
        # c = f * (h + e - 1)
        # Heights are additive in canonical towers
        # So h0 = c0 / f0 + 1 - e0
        if heights:
            self.c = f * (heights[-1] + e - 1)
            self.c_absolute = f * (heights[-1] + c0 / f0 - e0 + e)
        else:
            self.c = n - f
            self.c_absolute = f * (c0 / f0 - e0 + e)

    @lazy_attribute
    def scaled_heights(self):
        p, etame = self.p, self.etame
        return [h / (etame * p**i) for (i, h) in enumerate(self.heights, 1)]

    @lazy_attribute
    def scaled_rams(self):
        p, etame = self.p, self.etame
        return [r / (etame * p**i) for (i, r) in enumerate(self.rams, 1)]

    @lazy_attribute
    def top_slope(self):
        if self.slopes:
            return self.slopes[-1]
        elif self.e > 1:
            return 0
        else:
            return -1

    @lazy_attribute
    def slope_multiplicities(self):
        if not self.slopes:
            return []
        mult = [1]
        for s, t in zip(self.slopes[:-1], self.slopes[1:]):
            if s == t:
                mult[-1] += 1
            else:
                mult.append(1)
        return mult

    @lazy_attribute
    def wild_segments(self):
        return len(set(self.slopes))

    @lazy_attribute
    def bands(self):
        return [((0, 1+h), (self.e, h), (0, 1+s), (self.e, s)) for (h, s) in zip(self.scaled_heights, self.slopes)]

    @lazy_attribute
    def black(self):
        return [(0, 1), (self.e, 0)]

    @lazy_attribute
    def virtual_green(self):
        p, e, w = self.p, self.e, self.w
        last_slope = {}
        for i, s in enumerate(self.slopes, 1):
            last_slope[s] = i
        ans = []
        for i, (h, s) in enumerate(zip(self.scaled_heights, self.slopes), 1):
            u = e*frac(h)
            v = 1 + floor(h)
            if last_slope[s] == i:
                if (e*frac(h)).valuation(p) == w - i: # might be wrong
                    code = 1
                else:
                    code = 0
            else:
                code = -1
            ans.append((u, v, code))
        return ans

    @lazy_attribute
    def green(self):
        return [(u, v, bool(code)) for (u, v, code) in self.virtual_green if code >= 0]

    @lazy_attribute
    def solid_green(self):
        return [(u, v) for (u, v, solid) in self.green if solid]

    def _set_redblue(self):
        self.blue = []
        self.red = []
        p, e, w = self.p, self.e, self.w
        for i, (s, (u, v, code)) in enumerate(zip(self.slopes, self.virtual_green), 1):
            if u.denominator() == 1 and code == -1:
                self.blue.append((u, v, True))
            u = floor(u + 1)
            #print("Starting", i, s, u, v, code, e, w)
            while v <= 1 + s - u/e:
                #print("While", u, v, 1 + s - u/e, u.valuation(p), w-i)
                if u == e:
                    u = ZZ(0)
                    v += 1
                if v == 1 + s - u/e:
                    self.red.append((u, v, False))
                elif u.valuation(p) == (w - i):
                    self.blue.append((u, v, True))
                u += 1
        self.blue = sorted(set(self.blue))
        self.red = sorted(set(self.red))

    @lazy_attribute
    def blue(self):
        self._set_redblue()
        return self.blue

    @lazy_attribute
    def red(self):
        self._set_redblue()
        return self.red

    @lazy_attribute
    def poly(self):
        p, f = self.p, self.f
        pts = ([("a", u, v) for (u, v) in self.solid_green] +
               [("b", u, v) for (u, v, solid) in self.blue] +
               [("c", u, v) for (u, v, solid) in self.red])
        names = [f"{c}{self.e*(v-1)+u}" for (c, u, v) in pts]
        if gcd(p**f - 1, self.etame) > 1:
            names.append("d")
        if self.e0 > 1:
            names.append("pi")
        R = PolynomialRing(ZZ, names)
        if self.f == 1:
            S = PolynomialRing(R, "x")
        else:
            S = PolynomialRing(R, "nu")
        if self.e == 1:
            return S.gen()
        if "d" in names:
            d = R.gen(names.index("d"))
        else:
            d = 1
        x = S.gen()
        if self.e0 > 1:
            pi = R.gens()[-1]
        else:
            pi = p
        poly = x**self.e + d*pi
        for i, (c, u, v) in enumerate(pts):
            poly += R.gen(i) * pi**v * x**u
        return poly

    @lazy_attribute
    def ramification_polygon_old(self):
        # Old-style for now
        p = self.p
        L = [(self.e, 0)]
        if self.pw != self.e:
            L.append((self.pw, 0))
        cur = (self.pw, 0)
        for r, nextr in zip(self.rams, self.rams[1:] + [None]):
            x = cur[0] // p
            y = cur[1] + x * r
            cur = (x, y)
            if r != nextr:
                L.append(cur)
        L.reverse()
        return L

    @lazy_attribute
    def ramification_polygon_new(self):
        indices = self.generic_indices_insep
        p, w, pw, e, etame = self.p, self.w, self.pw, self.e, self.etame
        assert len(indices) == w + 1
        corners = [(-pw, 0)]
        edges = []
        a = w
        while a > 0:
            smallslope = min((indices[i] - indices[a]) / (p**a - p**i) for i in range(a))
            newedge = [corners[-1]]
            for i in range(1, a + 1):
                if (indices[a - i] - indices[a]) / (p**a - p**(a-i)) == smallslope:
                    newedge.append((-p**(a-i), indices[a-i]))
            edges.append(newedge)
            corners.append(newedge[-1])
            a = corners[-1][0].valuation(p)
        if etame > 1:
            corners.insert(0, (-e, 0))
            edges.insert(0, [(-e, 0), (-pw, 0)])
        return corners, edges

    @lazy_attribute
    def residual_polynomials(self):
        assert self.n0 == 1
        p, e = self.p, self.e
        if e == 1:
            return []
        f = self.poly
        S = f.parent()
        x = S.gen()
        Rp = S.base_ring().change_ring(GF(p)) # should maybe be GF(q)
        Sp = Rp['z']
        z = Sp.gen()

        def mpoly_valuation(g):
            pval, loc = min((min(c.valuation(p) for c in gi.coefficients()), i) for (i, gi) in g.dict().items())
            return e * pval + loc

        def pRed(g):
            if g == 0:
                return g
            pval, loc = min((min(c.valuation(p) for c in gi.coefficients()), i) for (i, gi) in g.dict().items())
            pred = p**(pval + 1)
            return S({i: {mon: c % pred for (mon, c) in gi.dict().items()} for (i, gi) in g.dict().items() if i <= loc})

        @cached_function
        def ramcoeffs(i):
            return pRed(S([j.binomial(i) * f[j] for j in srange(e + 1)]) % f)

        respolys = []
        verts = self.ramification_polygon_old
        for (x0, y0), (x1, y1) in zip(reversed(verts[:-1]), reversed(verts[1:])):
            xdiff = x1 - x0
            slope = (y0 - y1) / xdiff # notice negation
            a, b = slope.numerator(), slope.denominator()
            valx0 = mpoly_valuation(ramcoeffs(x0))
            newrespoly = Sp(0)
            for j in range((xdiff // b) + 1):
                ramcoeff = ramcoeffs(j*b + x0)
                xpower = pRed((x**(valx0 - j*a)) % f)
                vr = mpoly_valuation(ramcoeff)
                vx = mpoly_valuation(xpower)
                assert vx <= vr, ("Nonintegral quotient", ramcoeff, xpower)
                if vx == vr:
                    dr = ramcoeff.degree()
                    xcoeff, rcoeff = ZZ(xpower.coefficients()[0].coefficients()[0]), ZZ(ramcoeff.coefficients()[0].coefficients()[0])
                    frac = (rcoeff / xcoeff) % p
                    quot = frac * (ramcoeff // (rcoeff * x**dr))
                    newrespoly += quot[0].change_ring(GF(p)) * z**j
            respolys.append(newrespoly)
        return respolys

    @lazy_attribute
    def generic_indices_insep(self):
        """
        The indices of inseparability for the extension defined by substituting 1 for all the generic variables in this family
        """
        p, e, e0, w = self.p, self.e, self.e0, self.w
        f_generic = self.poly
        gens = f_generic.base_ring().gens()
        if e0 > 1:
            # need to keep pi
            Zx = ZZ['pi']['x']
            gens = gens[:-1]
            def valuation(c):
                if c == 0:
                    return infinity
                return min(i + e0*d.valuation(p) for (i, d) in c.dict())
        else:
            Zx = ZZ['x']
            def valuation(c):
                return c.valuation(p)
        D = {g: 1 for g in gens}
        f = Zx({i: c.subs(D) for (i, c) in f_generic.dict().items()})
        vals = [valuation(c) for c in f]
        tilde = [0] * (w + 1)
        indices = [0] * (w + 1)
        for j in range(w):
            tilde[w - j - 1] = min(i + e*vals[i] - e for i in srange(1, e) if i.valuation(p) <= w - j - 1)
            indices[w - j - 1] = min(tilde[w - j - 1], indices[w - j] + e*e0)
        return indices

    @lazy_attribute
    def residual_polynomials2(self):
        p, f, f0, e, etame = self.p, self.f, self.f0, self.e, self.etame
        _, edges = self.ramification_polygon_new
        if etame > 1:
            edges = edges[1:]
        poly = self.poly
        S = poly.parent()
        R = S.base_ring()
        names = tuple(poly.base_ring().variable_names())
        if "pi" in names:
            pi = R.gens()[-1]
            names = names[:-1]
        else:
            pi = p
        if "d" in names:
            names = names[:-1] + ("dinv",)
        q = p**(f * f0)
        Rp = GF(q)[names]
        Sp = Rp['z']
        z = Sp.gen()
        if "dinv" in names:
            mdinv = -Rp.gens()[-1]
        else:
            mdinv = -Rp(poly[0] // pi)
        respolys = []
        for E in edges:
            rpoly = {}
            p0 = -E[-1][0]
            slope = (E[-1][1] - E[0][1]) / (E[-1][0] - E[0][0])
            den = slope.denominator()
            for v in E:
                pj = -v[0]
                aj = (v[1] - 1) // e
                bj = v[1] - e * aj
                assert (pj - p0) % den == 0
                rpoly[(pj-p0)//den] = Rp((binomial(bj, pj) * poly[bj]) // pi**(aj + 1)) * mdinv**(aj + 1)
            respolys.append(Sp(rpoly))
        if etame > 1:
            tameres = (z + 1)**e - 1
            tameres //= z**(tameres.valuation())
            respolys.insert(0, tameres)
        return respolys

    @lazy_attribute
    def possible_phi0(self):
        assert self.f0 == self.f == 1 # Makes life simpler
        p, e0 = self.p, self.e0
        R = Zmod(p)
        possreps = set(range(1,p))
        e0pows = set(R(a)**e0 for a in possreps)
        minreps = []
        for a in range(1,p):
            if a in possreps:
                minreps.append(a)

    #@lazy_attribute
    #def possible_residual_polynomials(self):
    #    n = self.n
    #    for phi0 in 

    @lazy_attribute
    def gamma(self):
        if self.f_absolute == 1:
            return len(self.red)
        e = self.e
        gamma = 0
        for (u, v, _) in self.red:
            s = v + u/e - 1
            cnt = self.slopes.count(s)
            gamma += gcd(cnt, self.f_absolute)
        return gamma

    @lazy_attribute
    def poly_count(self):
        p, f, alpha, beta, gamma = self.p, self.f_absolute, len(self.solid_green), len(self.blue), self.gamma
        q = p**f
        return (q-1)**alpha * q**beta * p**gamma

    @lazy_attribute
    def mass_stored(self):
        return sum([ZZ(1) / rec["aut"] for rec in self.fields], ZZ(0))

    @lazy_attribute
    def mass(self):
        p, e0, f, alpha, beta = self.p, self.e0, self.f_absolute, len(self.solid_green), len(self.blue)
        q = p**f
        return (q-1)**alpha * q**beta / (self.f * self.base_aut)

    @lazy_attribute
    def mass_missing(self):
        return float((self.mass - self.mass_stored) / self.mass)

    @lazy_attribute
    def mass_display(self):
        return self.mass
        #m = self.mass
        #if m.denominator() == 1:
        #    return m
        #if m.floor() == 0:
        #    return latex(m)
        #return f"{m.floor()} {latex(m-m.floor())}"

    @lazy_attribute
    def ambiguity(self):
        return self.p**self.gamma

    @lazy_attribute
    def all_stored(self):
        return "t" if self.mass == self.mass_stored else "f"

    def respoly_iter(self):
        assert self.n0 == 1
        rp = self.residual_polynomials
        if rp:
            Rp = rp[0].base_ring()
            k = Rp.base_ring()
            opts = {}
            for f in rp:
                for c in f.coefficients():
                    for v in c.variables():
                        if v not in opts:
                            s = str(v)
                            if s[0] == "a":
                                opts[v] = [a for a in k if a != 0]
                            elif s[0] == "b":
                                opts[v] = list(k)
                            else:
                                raise ValueError("c found", self.rams, rp)
            kz = k['z']
            D = {v: k(0) for v in Rp.gens()}
            for vec in cartesian_product([[(v, a) for a in opts[v]] for v in opts]):
                D.update(dict(vec))
                yield [kz([c.subs(D) for c in f]) for f in rp]
        else:
            yield []

    def respoly2_iter(self):
        rp = self.residual_polynomials2
        seen = set()
        if rp:
            Rp = rp[0].base_ring()
            k = Rp.base_ring()
            q = k.cardinality()
            u = gcd(q - 1, self.etame)
            if u > 1:
                if q == self.p:
                    # Totally ramified
                    for a in range(2,q):
                        a = k(a)
                        if a.multiplicative_order() == q - 1:
                            break
                else:
                    a = k.gen()
                    assert a.multiplicative_order() == q - 1
                dinvs = [a**-i for i in range(u)]
            opts = {}
            for rpoly in rp:
                for c in rpoly.coefficients():
                    for v in c.variables():
                        if v not in opts:
                            s = str(v)
                            if s[0] == "a":
                                opts[v] = [a for a in k if a != 0]
                            elif s[0] == "b":
                                opts[v] = list(k)
                            elif s == "dinv":
                                opts[v] = dinvs
                            else:
                                raise ValueError("c found", self.rams, rp)
            kz = k['z']
            D = {v: k(0) for v in Rp.gens()}
            for vec in cartesian_product([[(v, a) for a in opts[v]] for v in opts]):
                D.update(dict(vec))
                new = [kz([c.subs(D) for c in rpoly]) for rpoly in rp]
                ntup = tuple(new)
                if ntup not in seen:
                    seen.add(ntup)
                    yield new
        else:
            yield []

    def respoly2_fixed_iter(self):
        # Convert to current residual polynomial convention
        _, edges = self.ramification_polygon_new
        q = self.p ** (self.f * self.f0)
        kz = GF(q)['z']
        for rp in self.respoly2_iter():
            fixed = []
            assert len(edges) == len(rp)
            for E, rpoly in zip(edges, rp):
                start = -E[-1][0]
                slope = (E[1][1] - E[0][1]) / (E[1][0] - E[0][0])
                a, b = slope.numerator(), slope.denominator()
                D = rpoly.dict()
                assert all((i - start) % b == 0 for i in D)
                fixed.append(kz({(i - start) // b: c for (i, c) in D.items()}))
            yield fixed

    def __iter__(self):
        assert self.n0 == 1
        f, e, etame, generic = self.f, self.e, self.etame, self.poly
        S = generic.parent()
        R = S.base_ring()
        Zx = PolynomialRing(ZZ, "x")
        x = Zx.gen()
        names = R.variable_names()
        p = self.p
        q = p**f
        opts = {}
        if f > 1:
            nu = self.nu
            k = GF(q, 'x', nu)
            d0 = k.gen()
            bopts = [Zx(y.polynomial()) for y in k]
            # We should really find the generic residual polynomial, specialize at a given
            # set of a and b values, compute Hermite normal form and pick representatives
            # for the quotient to find copts.  But this is annoying, so instead
            # we accept duplication since the bulk of the cases will have f=1 anyway
            #for u, v, _ in self.red:
            #    s = v + u/e - 1
            #    cnt = self.slopes.count(s)
        else:
            bopts = [ZZ(b) for b in range(p)]
        assert bopts[0] == 0
        aopts = bopts[1:]
        copts = bopts

        for name in names:
            if name[0] == "a":
                opts[name] = aopts
            elif name[0] == "b":
                opts[name] = bopts
            elif name[0] == "c":
                opts[name] = copts
        if self.w > 0 or "d" in names:
            if "d" in names:
                # Pick minimal d in the etame-th power classes in F_q.  When f>1 need to use nu...
                q = p**f
                dtop = gcd(q - 1, self.etame)
                if f > 1:
                    opts["d"] = [(d0**i).polynomial().change_ring(ZZ) for i in range(dtop)]
                else:
                    # q != 2 since dtop > 1
                    d0 = mod(2,p)
                    while d0.multiplicative_order() != p - 1:
                        d0 += 1
                    opts["d"] = [ZZ(d0**i) for i in range(dtop)]
            for vec in cartesian_product_iterator([opts[name] for name in names]):
                base_map = R.hom(vec)
                if f > 1:
                    Shom = S.hom([nu], base_map=base_map)
                else:
                    Shom = S.hom([x], base_map=base_map)
                yield Shom(generic)
        elif e > 1:
            # One possible tame extension
            return nu**e + p
        else:
            yield nu

    @lazy_attribute
    def cache_key(self):
        return self.p, self.f_absolute, self.e_absolute, str(self.visible)

    @lazy_attribute
    def field_query(self):
        return {"p": self.p,
                "f": self.f_absolute,
                "e": self.e_absolute,
                "visible": str(self.visible)}

    @lazy_attribute
    def fields(self):
        if self.field_cache is not None:
            L = self.field_cache[self.cache_key]
        else:
            L = list(db.lf_fields.search(self.field_query, cache_cols))
        if self.n0 == 1:
            return L
        # Prior version that only suppoted bases in the canonical filtration
        #return [rec for rec in L if self.base in rec["canonical_filtration"]]
        return [rec for rec in L if self.oldbase in rec["subfield"]]

    @lazy_attribute
    def nu(self):
        p, f = self.p, self.f
        q = p**f
        Zx = PolynomialRing(ZZ, "x")
        # Need to do something else if no conway polynomial is available
        k = GF(q, 'x', conway_polynomial(p, f))
        return Zx(k.polynomial())

    @lazy_attribute
    def field_sort_key(self):
        assert self.n0 == 1
        p, f, e, nu = self.p, self.f, self.e, self.nu
        def sort_key(poly):
            # poly should be in Eisenstein form with respect to self.nu
            if e == 1:
                # No sorting required
                return poly
            if f > 1:
                vec = []
                while poly.degree() > 0:
                    poly, r = poly.quo_rem(nu)
                    assert r % p == 0
                    vec.extend([r[j] // p for j in range(f)])
            else:
                assert all(c % p == 0 for c in list(poly)[:-1])
                vec = [c // p for c in list(poly)[:-1]]
            krasner_bound = (e * self.slopes[-1]).floor()
            i = 0
            ans = [0 for _ in range(f*(krasner_bound + 1))]
            while i <= krasner_bound:
                for j in range(f):
                    sig = (f*i + j)
                    bar = sig % len(vec)
                    vec[bar], ans[sig] = vec[bar].quo_rem(p)
                i += 1
            return ans
        return sort_key

    @lazy_attribute
    def field_count(self):
        return len(self.fields)

    @lazy_attribute
    def packet_count(self):
        if self.mass != self.mass_stored:
            return r"\N"
        return len(set(rec["packet"] for rec in self.fields))

    @classmethod
    def families(cls, base, n, f=None, field_cache=None, field_data=None):
        if isinstance(base, (Integer, int)):
            p, f0, e0 = base, ZZ(1), ZZ(1)
            base = f"{base}.1.1.0a1.1"
        else:
            p, f0, e0, _, _ = base.split(".")
            p, f0, e0 = ZZ(p), ZZ(f0), ZZ(e0)
        n0 = e0 * f0
        n = ZZ(n)
        if field_data is None:
            field_data = get_field_data()
        if field_cache is None:
            field_cache = defaultdict(list)
            for rec in db.lf_fields.search({"p": p, "n": n*n0}, cache_cols):
                #if n0 == 1 or base in rec["canonical_filtration"]:
                if n0 == 1 or field_data[base]["label"] in rec["subfield"]:
                    field_cache[cache_key(rec)].append(rec)
        if f is None:
            for f in reversed(n.divisors()):
                # maybe this should have base_aut = f?
                yield from cls.families(base, n, f, field_cache=field_cache, field_data=field_data)
        else:
            e = n // f
            w, etame = e.val_unit(p)
            if n0 == 1:
                base_data = {"aut": 1, "visible": [], "rf": [1,0], "rams": [], "label": f"{p}.1.0.1"}
            else:
                base_data = field_data[base]
                if base_data["visible"]:
                    w0, e0tame = e0.val_unit(p)
                    base_slopes = [s - 1 for s in base_data["visible"]]
                    base_data["rams"] = heights_to_rams(slopes_to_heights(base_slopes, p, e0tame), p)
                else:
                    base_data["rams"] = []

            def R(e, rho):
                #den = sum(p**j for j in range(rho))
                #nums = [n for n in range(1, p*e*den) if n % p != 0 and n > base_lastram * den]
                #if rho == 1 and base_lastram < p*e:
                #    nums.append(p*e*den)
                den = p**rho - 1
                nums = [n for n in range(1, p*e*den//(p-1)) if n % p != 0]
                if rho == 1:
                    nums.append(p*e*den//(p-1))
                return [n / den for n in nums]

            # We can't compute the label intrinsically, so we have to collect and sort
            fams = []
            if w == 0:
                fams.append(cls(base, f, etame, base_data=base_data, rams=[], field_cache=field_cache))
            else:
                for mvec in reversed(OrderedPartitions(w)):
                    Mvec = [0]
                    for m in mvec[:-1]:
                        Mvec.append(Mvec[-1] + m)
                    Rs = [R(e0 * etame * p**M, m) for m, M in zip(mvec, Mvec)]
                    for rvec in cartesian_product(Rs):
                        if all(a < b for (a,b) in zip(rvec[:-1], rvec[1:])):
                            rmvec = []
                            for r, m in zip(rvec, mvec):
                                rmvec.extend([r] * m)
                            fams.append(cls(base, f, etame, base_data=field_data[base], rams=rmvec, field_cache=field_cache))
                fams.sort(key=lambda fam: (fam.c, len(fam.slope_multiplicities), list(reversed(fam.slope_multiplicities)), fam.rams))
            ctr = Counter()
            for fam in fams:
                fam.ctr = ctr[fam.c]
                code = cremona_letter_code(fam.ctr)
                ctr[fam.c] += 1
                if n0 == 1:
                    # Absolute labels
                    fam.label = f"{p}.{f}.{e}.{fam.c}{code}"
                else:
                    # Relative labels
                    fam.label = f"{base}-{f}.{e}.{fam.c}{code}"
            yield from fams


