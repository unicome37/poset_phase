import json
from pathlib import Path

data = json.loads(Path('d:/Kiro/理论体系/poset_phase/outputs_carlip/f2_turnon_margin_refit.json').read_text(encoding='utf-8'))
os_ = data['onset_summary']
cn = os_.get('consecutive_n', 3)
for level in ['winner_only', 'operational', 'manuscript_safe']:
    info = os_.get(level, {})
    fbs = info.get('first_block_start', 'N/A')
    blk = info.get('blocks', [])
    print(f"{level}: first_block_start={fbs}, consecutive_n={cn}, block={blk[0] if blk else None}")
