import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =========================
# Leer trayectoria
# =========================
df = pd.read_csv("position.csv")
traj = df[["x", "y", "z"]].values
time = df["time"].values

# =========================
# Obstáculos reales
# =========================
obstacles = [
    {"name": "Obs 1", "center": np.array([1.1, -0.5, 1.0]), "radius": 0.8},
    {"name": "Obs 2", "center": np.array([2.6,  0.5, 1.1]), "radius": 0.8},
    {"name": "Obs 3", "center": np.array([4.0,  1.6, 1.2]), "radius": 0.8},
]

# =========================
# Distancia punto-segmento
# =========================
def point_to_segment_distance(center, p1, p2):
    seg = p2 - p1
    seg_norm_sq = np.dot(seg, seg)

    if seg_norm_sq == 0:
        closest = p1
    else:
        t = np.dot(center - p1, seg) / seg_norm_sq
        t = np.clip(t, 0.0, 1.0)
        closest = p1 + t * seg

    dist = np.linalg.norm(center - closest)
    return dist, closest

# =========================
# Verificación
# =========================
print("\n=== Verificación de colisión ===")

for obs in obstacles:
    name = obs["name"]
    center = obs["center"]
    radius = obs["radius"]

    # 1) revisión punto a punto
    point_dists = np.linalg.norm(traj - center, axis=1)
    min_point_dist = point_dists.min()
    idx_point = point_dists.argmin()
    clearance_points = min_point_dist - radius

    # 2) revisión por segmentos
    min_seg_dist = float("inf")
    closest_point = None
    closest_seg_idx = None

    for i in range(len(traj) - 1):
        d, cp = point_to_segment_distance(center, traj[i], traj[i + 1])
        if d < min_seg_dist:
            min_seg_dist = d
            closest_point = cp
            closest_seg_idx = i

    clearance_seg = min_seg_dist - radius

    print(f"\n{name}")
    print(f"  Distancia mínima por puntos   : {min_point_dist:.6f}")
    print(f"  Clearance por puntos          : {clearance_points:.6f}")
    print(f"  Tiempo aprox punto más cercano: {time[idx_point]:.2f} s")
    print(f"  Punto más cercano             : {traj[idx_point]}")

    print(f"  Distancia mínima por segmento : {min_seg_dist:.6f}")
    print(f"  Clearance por segmento        : {clearance_seg:.6f}")
    print(f"  Segmento más cercano          : entre t={time[closest_seg_idx]:.2f}s y t={time[closest_seg_idx+1]:.2f}s")
    print(f"  Punto más cercano en segmento : {closest_point}")

    if min_seg_dist < radius:
        print("  >>> COLISIÓN: la trayectoria entra en el obstáculo")
    elif np.isclose(min_seg_dist, radius, atol=1e-6):
        print("  >>> ROZA el obstáculo")
    else:
        print("  >>> No colisiona")

# =========================
# Gráfica 3D
# =========================
fig = plt.figure(figsize=(9, 7))
ax = fig.add_subplot(111, projection="3d")

x = traj[:, 0]
y = traj[:, 1]
z = traj[:, 2]

ax.plot(x, y, z, linewidth=2.5, label="Trayectoria")
ax.scatter(x[0], y[0], z[0], s=70, label="Inicio")
ax.scatter(x[-1], y[-1], z[-1], s=70, label="Final")

u = np.linspace(0, 2*np.pi, 40)
v = np.linspace(0, np.pi, 25)

for obs in obstacles:
    cx, cy, cz = obs["center"]
    r = obs["radius"]

    xs = cx + r * np.outer(np.cos(u), np.sin(v))
    ys = cy + r * np.outer(np.sin(u), np.sin(v))
    zs = cz + r * np.outer(np.ones_like(u), np.cos(v))

    ax.plot_surface(xs, ys, zs, alpha=0.25)
    ax.text(cx, cy, cz, obs["name"])

ax.set_xlabel("X [m]")
ax.set_ylabel("Y [m]")
ax.set_zlabel("Z [m]")
ax.set_title("Trayectoria 3D con obstáculos")
ax.legend()

max_range = max(
    x.max() - x.min(),
    y.max() - y.min(),
    z.max() - z.min()
)

mid_x = (x.max() + x.min()) / 2
mid_y = (y.max() + y.min()) / 2
mid_z = (z.max() + z.min()) / 2

ax.set_xlim(mid_x - max_range/2, mid_x + max_range/2)
ax.set_ylim(mid_y - max_range/2, mid_y + max_range/2)
ax.set_zlim(mid_z - max_range/2, mid_z + max_range/2)

plt.tight_layout()
plt.show()