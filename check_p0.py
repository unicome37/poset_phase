import json
from pathlib import Path

base = Path('d:/Kiro/理论体系/poset_phase/outputs_carlip')

for name in ['falsify_c2_background_response_flrw_metricfaithful_phaseA',
             'falsify_c2_background_response_flrw_proxy_control_phaseA']:
    partial = base / f'{name}.partial.json'
    final = base / f'{name}.json'
    
    if final.exists():
        data = json.loads(final.read_text(encoding='utf-8'))
        status = 'FINAL'
    elif partial.exists():
        data = json.loads(partial.read_text(encoding='utf-8'))
        status = 'PARTIAL'
    else:
        print(name, '=> NOT FOUND')
        continue
    
    label = 'METRIC' if 'metric' in name else 'PROXY'
    recs = data.get('records', [])
    done = data.get('completed_seed_runs', 0)
    total = data.get('total_seed_runs', '?')
    print(f'=== {label} [{status}]: {done}/{total} seeds, {len(recs)} records ===')
    for r in recs:
        seed_r = r.get('seed_run', '?')
        N_r = r.get('N', '?')
        param = r.get('param', '?')
        rank = r.get('rank', '?')
        winner = r.get('winner', 'N/A')
        print(f'  seed={seed_r} N={N_r} kappa={param} rank={rank} winner={winner}')
    print()
