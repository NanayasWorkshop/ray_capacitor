# capacitor_viewer.py
import numpy as np
import trimesh
import sys
import os

def load_data(filename):
    """Load ray tracing data from exported text file"""
    if not os.path.exists(filename):
        filename = os.path.join('..', filename)
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    pos_vertices = []
    pos_faces = []
    neg_vertices = []
    neg_faces = []
    rays = []
    
    mode = None
    for line in lines:
        line = line.strip()
        if line == 'POSITIVE_MESH':
            mode = 'pos'
        elif line == 'NEGATIVE_MESH':
            mode = 'neg'
        elif line == 'RAYS':
            mode = 'rays'
        elif line:
            vals = list(map(float, line.split()))
            if mode == 'pos':
                v_idx = len(pos_vertices)
                pos_vertices.extend([vals[0:3], vals[3:6], vals[6:9]])
                pos_faces.append([v_idx, v_idx+1, v_idx+2])
            elif mode == 'neg':
                v_idx = len(neg_vertices)
                neg_vertices.extend([vals[0:3], vals[3:6], vals[6:9]])
                neg_faces.append([v_idx, v_idx+1, v_idx+2])
            elif mode == 'rays':
                rays.append(vals)
    
    return np.array(pos_vertices), np.array(pos_faces), \
           np.array(neg_vertices), np.array(neg_faces), rays

def create_visualization(pos_verts, pos_faces, neg_verts, neg_faces, rays):
    """Create trimesh visualization"""
    
    # Determine which faces emit rays (for positive mesh)
    ray_origins = np.array([ray[0:3] for ray in rays])
    
    # Calculate face centers for positive mesh
    pos_face_centers = pos_verts[pos_faces].mean(axis=1)
    
    # Color positive mesh - darker red for emitting faces
    pos_colors = np.zeros((len(pos_faces), 4))
    for i, center in enumerate(pos_face_centers):
        # Check if this face is near any ray origin
        distances = np.linalg.norm(ray_origins - center, axis=1)
        if np.any(distances < 0.5):  # Emitting face
            pos_colors[i] = [128, 0, 0, 200]  # Dark red
        else:
            pos_colors[i] = [255, 100, 100, 150]  # Light red
    
    # Create positive mesh
    pos_mesh = trimesh.Trimesh(vertices=pos_verts, faces=pos_faces)
    pos_mesh.visual.face_colors = pos_colors
    
    # Calculate which negative faces are hit
    hit_endpoints = np.array([ray[3:6] for ray in rays if ray[6] > 0])
    neg_face_centers = neg_verts[neg_faces].mean(axis=1)
    
    # Color negative mesh - darker blue for hit faces
    neg_colors = np.zeros((len(neg_faces), 4))
    for i, center in enumerate(neg_face_centers):
        if len(hit_endpoints) > 0:
            distances = np.linalg.norm(hit_endpoints[:, None] - center, axis=2).min()
            if distances < 0.5:  # Hit face
                neg_colors[i] = [0, 0, 128, 200]  # Dark blue
            else:
                neg_colors[i] = [100, 100, 255, 150]  # Light blue
        else:
            neg_colors[i] = [100, 100, 255, 150]  # Light blue
    
    # Create negative mesh
    neg_mesh = trimesh.Trimesh(vertices=neg_verts, faces=neg_faces)
    neg_mesh.visual.face_colors = neg_colors
    
    # Create ray visualizations - only every 10th ray
    ray_paths = []
    hits = 0
    misses = 0
    
    print(f"Processing {len(rays)} rays, displaying every 10th...")
    
    for i, ray in enumerate(rays):
        is_hit = ray[6] > 0
        
        if is_hit:
            hits += 1
        else:
            misses += 1
        
        # Only show every 10th ray
        if i % 10 != 0:
            continue
            
        origin = ray[0:3]
        endpoint = ray[3:6]
        
        if is_hit:
            # Create green line for hits
            path = trimesh.load_path([origin, endpoint])
            path.colors = [[0, 255, 0, 255]]  # Green
            ray_paths.append(path)
        else:
            # Create short orange line for misses
            direction = np.array(endpoint) - np.array(origin)
            norm = np.linalg.norm(direction)
            if norm > 0:
                direction = direction / norm * 2.0  # 2mm for visualization
            endpoint_short = np.array(origin) + direction
            path = trimesh.load_path([origin, endpoint_short])
            path.colors = [[255, 165, 0, 180]]  # Orange
            ray_paths.append(path)
    
    # Combine all geometry into scene
    scene = trimesh.Scene()
    scene.add_geometry(pos_mesh, node_name='positive_plate')
    scene.add_geometry(neg_mesh, node_name='negative_plate')
    
    # Add rays to scene
    for i, path in enumerate(ray_paths):
        scene.add_geometry(path, node_name=f'ray_{i}')
    
    print(f"\nVisualization Statistics:")
    print(f"Positive mesh: {len(pos_faces)} triangles")
    print(f"  - Emitting faces (dark red): ray origins")
    print(f"  - Non-emitting faces (light red): other faces")
    print(f"Negative mesh: {len(neg_faces)} triangles")
    print(f"  - Target faces (dark blue): ray endpoints")
    print(f"  - Non-target faces (light blue): other faces")
    print(f"Rays: {hits} hits (green), {misses} misses (orange)")
    print(f"Showing: {len(ray_paths)} rays (every 10th)")
    print(f"\nControls:")
    print("- Mouse drag: Rotate view")
    print("- Scroll: Zoom")
    print("- W: Toggle wireframe")
    print("- C: Toggle backface culling")
    print("- A: Toggle axis")
    print("- F: Reset view")
    
    # Show the scene
    scene.show()

def main():
    print("Capacitor Ray Tracing Visualizer (Trimesh)")
    print("-" * 40)
    
    try:
        pos_verts, pos_faces, neg_verts, neg_faces, rays = load_data('ray_data.txt')
        print(f"Loaded data successfully")
    except FileNotFoundError:
        print("Error: ray_data.txt not found!")
        print("Make sure to run: ray_distance.exe -v")
        sys.exit(1)
    except Exception as e:
        print(f"Error loading data: {e}")
        sys.exit(1)
    
    create_visualization(pos_verts, pos_faces, neg_verts, neg_faces, rays)

if __name__ == "__main__":
    main()