import numpy as np
import matplotlib.pyplot as plt

# Ler o arquivo points.obj
points = []
with open('points.obj', 'r') as f:
    for line in f:
        if line.startswith('v '):
            parts = line.strip().split()
            if len(parts) >= 4:
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])
                points.append([x, y, z])

points = np.array(points)
print(f"Total de pontos: {len(points):,}")

# ==================== SCATTER PLOT SIMPLES ====================
plt.figure(figsize=(10, 9))

# Projeção XZ (melhor vista da espiral)
plt.scatter(points[:, 0], points[:, 2], 
            c=points[:, 1], 
            cmap='viridis', 
            s=3, 
            alpha=0.8)

plt.title('CORDIC Spiral - Scatter Plot Simples (Projeção XZ)', fontsize=14)
plt.xlabel('X')
plt.ylabel('Z')
plt.axis('equal')
plt.grid(True, alpha=0.3)

# Barra de cores
plt.colorbar(label='Eixo Y')

plt.tight_layout()
plt.savefig('cordic_scatter_simples.png', dpi=300)
print("✅ Imagem salva como 'cordic_scatter_simples.png'")

# plt.show()   # descomente se quiser ver na tela