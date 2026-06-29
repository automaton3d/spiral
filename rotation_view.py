"""
rotation_view.py — Interactive 3D visualization of before/after rotation.

Usage:
  python rotation_view.py                       # default files
  python rotation_view.py before.dat after.dat   # custom files

Requires: matplotlib, numpy
  pip install matplotlib numpy
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

# Subsample for performance if too many points
MAX_PTS = 50000
if len(before) > MAX_PTS:
    idx = np.random.choice(len(before), MAX_PTS, replace=False)
    before_show = before[idx]
else:
    before_show = before

if len(after) > MAX_PTS:
    idx = np.random.choice(len(after), MAX_PTS, replace=False)
    after_show = after[idx]
else:
    after_show = after

fig = plt.figure(figsize=(16, 7))
fig.suptitle(f'Rotation  (L={L}, {len(before)} bits)', fontsize=14)

# Before — 3D scatter
ax1 = fig.add_subplot(121, projection='3d')
ax1.scatter(before_show[:, 0], before_show[:, 1], before_show[:, 2],
            s=0.3, c='blue', alpha=0.4, depthshade=True)
ax1.set_xlim(0, L); ax1.set_ylim(0, L); ax1.set_zlim(0, L)
ax1.set_xlabel('X'); ax1.set_ylabel('Y'); ax1.set_zlabel('Z')
ax1.set_title(f'Before ({len(before)} pts)')
ax1.set_box_aspect([1, 1, 1])

# After — 3D scatter
ax2 = fig.add_subplot(122, projection='3d')
ax2.scatter(after_show[:, 0], after_show[:, 1], after_show[:, 2],
            s=0.3, c='red', alpha=0.4, depthshade=True)
ax2.set_xlim(0, L); ax2.set_ylim(0, L); ax2.set_zlim(0, L)
ax2.set_xlabel('X'); ax2.set_ylabel('Y'); ax2.set_zlabel('Z')
ax2.set_title(f'After ({len(after)} pts)')
ax2.set_box_aspect([1, 1, 1])

# Sync rotation between the two panels
def on_move(event):
    if event.inaxes == ax1:
        ax2.view_init(elev=ax1.elev, azim=ax1.azim)
    elif event.inaxes == ax2:
        ax1.view_init(elev=ax2.elev, azim=ax2.azim)
    fig.canvas.draw_idle()

fig.canvas.mpl_connect('motion_notify_event', on_move)

plt.tight_layout()
print(f'Grid: {L}x{L}x{L}')
print(f'Bits before: {len(before)}, after: {len(after)}, delta: {len(after)-len(before)}')
print('Rotate with mouse. Close window to exit.')
plt.show()
