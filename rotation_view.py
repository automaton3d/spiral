"""
rotation_view.py — Visualize before/after rotation from ac_rotation.c output.

Usage:
  python rotation_view.py                       # default files
  python rotation_view.py before.dat after.dat   # custom files
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load_state(filename):
    data = []
    L = 128
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('# L ='):
                L = int(line.split('=')[1].strip())
            elif line.startswith('#') or line.startswith('x'):
                continue
            else:
                parts = line.strip().split()
                if len(parts) >= 3:
                    data.append([int(parts[0]), int(parts[1]), int(parts[2])])
    return np.array(data), L

before_file = sys.argv[1] if len(sys.argv) > 1 else 'before_rotation.dat'
after_file  = sys.argv[2] if len(sys.argv) > 2 else 'after_rotation.dat'

before, L = load_state(before_file)
after, _  = load_state(after_file)

z_mid = L // 2

b_slice = before[before[:, 2] == z_mid]
a_slice = after[after[:, 2] == z_mid]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

ax1.scatter(b_slice[:, 0], b_slice[:, 1], s=1, c='blue', marker='s')
ax1.set_xlim(0, L)
ax1.set_ylim(0, L)
ax1.set_aspect('equal')
ax1.set_title(f'Before (z={z_mid}, {len(b_slice)} pts)')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.grid(True, alpha=0.2)

ax2.scatter(a_slice[:, 0], a_slice[:, 1], s=1, c='red', marker='s')
ax2.set_xlim(0, L)
ax2.set_ylim(0, L)
ax2.set_aspect('equal')
ax2.set_title(f'After rotation (z={z_mid}, {len(a_slice)} pts)')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.grid(True, alpha=0.2)

plt.tight_layout()
plt.savefig('rotation_comparison.png', dpi=150, bbox_inches='tight')
print(f'Saved: rotation_comparison.png')
print(f'Grid: {L}x{L}x{L}')
print(f'Bits before: {len(before)}, after: {len(after)}, delta: {len(after)-len(before)}')
