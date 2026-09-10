"""Read-only Phase4 saved-band audit; stdlib only, no calculator or optimizer.

Usage: python -B diagnose_phase4_neb.py ATTEMPT [ATTEMPT ...]
Prints JSON to stdout; never creates evidence or changes scientific acceptance.
Math uses nonperiodic Cartesian coordinates in A, energies in eV, forces in eV/A.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re


def dot(a, b):
    if len(a) != len(b):
        raise ValueError('Vector dimensions differ')
    return sum(x * y for x, y in zip(a, b))


def norm(a):
    return math.sqrt(dot(a, a))


def sub(a, b):
    if len(a) != len(b):
        raise ValueError('Vector dimensions differ')
    return [x - y for x, y in zip(a, b)]


def finite(values):
    if not all(math.isfinite(x) for x in values):
        raise ValueError('Nonfinite input')


def atom_fmax(f):
    if not f or len(f) % 3:
        raise ValueError('Expected flattened Cartesian atom forces')
    finite(f)
    return max(norm(f[i:i + 3]) for i in range(0, len(f), 3))


def improved_tangent(backward, forward, energies):
    """ASE ImprovedTangentMethod convention, rejecting undefined directions."""
    if len(backward) != len(forward) or len(energies) != 3:
        raise ValueError('Invalid tangent dimensions')
    finite([*backward, *forward, *energies])
    left, center, right = energies
    if right > center > left:
        tangent = list(forward)
    elif right < center < left:
        tangent = list(backward)
    else:
        hi = max(abs(right - center), abs(left - center))
        lo = min(abs(right - center), abs(left - center))
        fw, bw = (hi, lo) if right > left else (lo, hi)
        tangent = [fw * f + bw * b for b, f in zip(backward, forward)]
    length = norm(tangent)
    if length <= 1e-14:
        raise ValueError('Undefined improved tangent (flat energies or degenerate path)')
    return [x / length for x in tangent]


def project_force(force, tangent, spring_parallel=0.0, climb=False):
    """Force must already have any ASE constraint adjustment applied."""
    finite([*force, *tangent, spring_parallel])
    if not math.isclose(norm(tangent), 1.0, abs_tol=1e-12):
        raise ValueError('Tangent must be unit normalized across ALL atoms')
    parallel = dot(force, tangent)
    perpendicular = [f - parallel * t for f, t in zip(force, tangent)]
    coefficient = -parallel if climb else spring_parallel
    projected = [f + coefficient * t for f, t in zip(perpendicular, tangent)]
    return {'projected': projected, 'perpendicular': perpendicular,
            'parallel_scalar': parallel, 'fmax': atom_fmax(projected),
            'physical_fmax': atom_fmax(force)}


def read_xyz(path):
    lines = path.read_text(encoding='utf-8').splitlines()
    frames, identity, cursor = [], None, 0
    while cursor < len(lines):
        if not lines[cursor].strip():
            cursor += 1
            continue
        n = int(lines[cursor])
        if n < 1 or cursor + n + 2 > len(lines):
            raise ValueError('Truncated/empty XYZ frame')
        symbols, positions = [], []
        for row in lines[cursor + 2:cursor + n + 2]:
            fields = row.split()
            if len(fields) != 4:
                raise ValueError('Only species + Cartesian position XYZ supported')
            symbols.append(fields[0])
            positions.extend(float(x) for x in fields[1:])
        finite(positions)
        if identity is not None and symbols != identity:
            raise ValueError('Atom order/identity differs between images')
        identity = symbols
        frames.append(positions)
        cursor += n + 2
    if not frames:
        raise ValueError('No XYZ frames')
    return identity, frames


def band_geometry(positions, energies, spring=0.1):
    if len(positions) < 3 or len(positions) != len(energies):
        raise ValueError('Need matching band and energies with interior images')
    size = len(positions[0])
    if not size or size % 3 or any(len(p) != size for p in positions):
        raise ValueError('Invalid band coordinate dimensions')
    finite([*energies, spring, *(v for p in positions for v in p)])
    if spring < 0:
        raise ValueError('Negative spring')
    highest = max(energies[1:-1])
    candidates = [i for i in range(1, len(energies) - 1) if energies[i] == highest]
    # Do not silently choose an ASE argsort-dependent tie.
    climbing = candidates[0] if len(candidates) == 1 else None
    segments = [sub(b, a) for a, b in zip(positions, positions[1:])]
    lengths = [norm(d) for d in segments]
    rows = []
    for i in range(1, len(positions) - 1):
        scalar = spring * (lengths[i] - lengths[i - 1])
        row = {'image': i, 'climbing': i == climbing,
               'spring_parallel_eV_A': scalar,
               'applied_spring_parallel_eV_A': 0.0 if i == climbing else scalar}
        try:
            tangent = improved_tangent(segments[i - 1], segments[i], energies[i - 1:i + 2])
            row['tangent'] = tangent
            # ||Fperp+s*t|| >= |s|; max atom norm >= band norm/sqrt(N).
            row['projected_fmax_lower_bound_eV_A'] = (
                0.0 if i in candidates else abs(scalar) / math.sqrt(size // 3))
        except ValueError as exc:
            row['tangent_error'] = str(exc)
        rows.append(row)
    contacts = []
    for p in positions:
        contacts.append(min(norm(sub(p[a:a + 3], p[b:b + 3]))
                            for a in range(0, size, 3) for b in range(a + 3, size, 3))
                        if size > 3 else None)
    return {'climbing_image': climbing, 'climbing_candidates': candidates,
            'segment_lengths_A': lengths, 'minimum_pair_distances_A': contacts,
            'images': rows}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replay_force_diagnostics(saved):
    """Reconstruct projected forces from same-state saved adjusted forces."""
    if saved['schema_version'] != 1 or saved['method'] != 'improvedtangent':
        raise ValueError('Unsupported force diagnostic schema/method')
    if saved['units'] != dict(positions='A', energies='eV', forces='eV/A', spring='eV/A^2'):
        raise ValueError('Unsupported force diagnostic units')
    if any(any(flags) for flags in saved['pbc']):
        raise ValueError('Periodic tangent replay is not implemented')
    flatten = lambda frame: [v for atom in frame for v in atom]
    positions = [flatten(p) for p in saved['positions_A']]
    energies = saved['energies_eV']
    geometry = band_geometry(positions, energies)
    if saved['force_image_indices'] != list(range(1, len(positions) - 1)):
        raise ValueError('Interior force indices do not match band')
    n = len(positions) - 2
    if any(len(saved[key]) != n for key in (
            'raw_forces_eV_A', 'constraint_adjusted_forces_eV_A', 'projected_forces_eV_A')):
        raise ValueError('Force image count mismatch')
    if len(saved['spring_eV_A2']) != len(positions) - 1:
        raise ValueError('Spring count mismatch')
    errors, maxima, constraint_changes = [], [], []
    for j, row in enumerate(geometry['images']):
        i = row['image']
        actual = flatten(saved['projected_forces_eV_A'][j])
        raw = flatten(saved['raw_forces_eV_A'][j])
        adjusted = flatten(saved['constraint_adjusted_forces_eV_A'][j])
        lengths, springs = geometry['segment_lengths_A'], saved['spring_eV_A2']
        result = project_force(adjusted, row['tangent'],
                               springs[i] * lengths[i] - springs[i - 1] * lengths[i - 1],
                               saved['climb'] and saved['climbing_image'] == i)
        errors.append(atom_fmax(sub(result['projected'], actual)))
        maxima.append(atom_fmax(actual))
        constraint_changes.append(atom_fmax(sub(adjusted, raw)))
    return dict(projection_max_atom_error_eV_A=max(errors),
                recomputed_saved_fmax_eV_A=max(maxima),
                fmax_bookkeeping_error_eV_A=max(maxima) - saved['projected_fmax_eV_A'],
                constraint_adjustment_fmax_eV_A=max(constraint_changes),
                force_gate_passes=max(maxima) <= saved['force_target_eV_A'],
                acceptance='diagnostic_only_no_TS_no_IRC')


def diagnose(attempt, spring=0.1):
    attempt = Path(attempt).resolve()
    out = attempt / 'results_phase4'
    report = {'attempt': str(attempt), 'acceptance': 'NOT_ESTABLISHED_NO_TS_NO_IRC',
              'units': {'coordinates': 'A', 'energy': 'eV', 'force': 'eV/A'},
              'spring_eV_A2': spring, 'spring_provenance': 'CLI assumption; verify run command',
              'sha256': {}, 'limitations': []}
    for path in [attempt / 'run.log', attempt / 'rerun_status.json',
                 attempt / 'run_phase4_reaction_mechanism.py',
                 *sorted(out.glob('*'))]:
        if path.is_file():
            report['sha256'][str(path.relative_to(attempt))] = sha256(path)
    if not (out / 'stage2.json').exists():
        report['status'] = 'incomplete_no_saved_final_band_metrics'
        return report
    record = json.loads((out / 'stage2.json').read_text())
    symbols, positions = read_xyz(out / 'neb_final_path.xyz')
    energies = record['energies_ev']
    geometry = band_geometry(positions, energies, spring)
    report['geometry'] = geometry
    report['atom_count'] = len(symbols)
    report['reported_fmax_eV_A'] = record['fmax_final_ev_A']
    target = record['fmax_target']
    finite([record['fmax_final_ev_A'], target])
    if target <= 0 or record['fmax_final_ev_A'] < 0:
        raise ValueError('Invalid force gate')
    force_pass = record['fmax_final_ev_A'] <= target
    report['target_eV_A'] = target
    report['bookkeeping'] = {
        'force_gate_passes': force_pass,
        'false_positive_convergence': bool(record['converged']) and not force_pass,
        'image_count_matches': record['n_images'] == len(positions),
        'ts_index_matches_interior_max': record['ts_image_index'] == geometry['climbing_image'],
        'forward_barrier_error_eV': record['barrier_fwd_ev'] - (max(energies) - energies[0])}
    initial_symbols, initial = read_xyz(out / 'images_idpp.xyz')
    report['endpoint_max_displacement_A'] = (
        max(atom_fmax(sub(initial[i], positions[i])) for i in (0, -1))
        if initial_symbols == symbols and len(initial) == len(positions) else None)
    log = (attempt / 'run.log').read_text(encoding='utf-8', errors='replace')
    report['export_results_attribute_error'] = (
        "'_PerAtomCalc' object has no attribute 'results'" in log)
    candidate = out / 'ts_candidate.xyz'
    report['candidate_bytes'] = candidate.stat().st_size if candidate.exists() else None
    report['logged_refine_fmax_eV_A'] = [float(x) for x in re.findall(r'fmax=([\d.]+) eV/A', log)]
    report['saved_force_files'] = [p.name for p in out.iterdir()
                                   if p.is_file() and ('force' in p.name.lower() or p.suffix == '.traj')]
    report['limitations'] = [
        'The inspected historical outputs contain positions/energies and scalar fmax, not raw or projected force arrays.',
        'Cannot independently reproduce reported fmax or assign residual to physical, tangent, projection, or constraint error.',
        'Geometry spacing and spring lower bounds are conditional diagnostics, not stationarity evidence.',
        'No per-step band history: cannot reconstruct the best-step index or explain FIRE divergence.',
        'Plain XYZ omits constraints/PBC/cell; inspected production source uses unconstrained nonperiodic molecular images.',
        'No QM, finite difference, Hessian, TS refinement, IRC, or NEB optimization was executed.']
    report['status'] = 'reported_neb_unconverged_force_cause_unresolved' if not force_pass else 'reported_force_gate_only'
    force_path = out / 'neb_force_diagnostics.json'
    if force_path.exists():
        saved = json.loads(force_path.read_text(encoding='utf-8'))
        report['force_replay'] = replay_force_diagnostics(saved)
        report['force_replay']['stage2_fmax_error_eV_A'] = (
            report['force_replay']['recomputed_saved_fmax_eV_A'] - record['fmax_final_ev_A'])
        report['force_replay']['stage2_energies_match'] = saved['energies_eV'] == energies
        saved_positions = [[v for atom in p for v in atom] for p in saved['positions_A']]
        report['force_replay']['xyz_matches_within_rounding'] = (
            len(saved_positions) == len(positions) and
            all(max(abs(v) for v in sub(a, b)) <= 5.1e-9
                for a, b in zip(saved_positions, positions)))
        report['limitations'] = report['limitations'][2:4] + [
            'Replay checks saved arrays, not independent QM forces or Hessian/IRC acceptance.',
            'No QM or NEB optimization was executed.']
        report['status'] = 'saved_force_replay_available_no_TS_no_IRC'
    # Keep tangents inspectable through band_geometry(), but shorten command output.
    for row in geometry['images']:
        if 'tangent' in row:
            row['tangent_norm'] = norm(row.pop('tangent'))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('attempts', nargs='+', type=Path)
    parser.add_argument('--spring', type=float, default=0.1,
                        help='Assumed scalar spring in eV/A^2; historical source default 0.1')
    args = parser.parse_args()
    print(json.dumps([diagnose(p, args.spring) for p in args.attempts], indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
