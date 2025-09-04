# multi_sensor_animated_viewer.py
import numpy as np
import trimesh
import sys
import os
import glob
import time
import math

class MultiSensorViewer:
    def __init__(self, radius=26.45):
        self.radius = radius
        self.sensor_positions = self.calculate_sensor_positions()
        self.time_steps = []
        self.current_step = 0
        self.scene = None
        
    def calculate_sensor_positions(self):
        """Calculate triangular positions for sensors A, B, C on ZX plane (around Y axis)"""
        positions = {}
        
        # Triangular formation with 120 degree spacing around Y axis (ZX plane)
        for i, sensor in enumerate(['A', 'C', 'B']):
            angle = i * (2 * math.pi / 3) + (math.pi / 2)  # Add 90° rotation
            x = self.radius * math.cos(angle)  # X coordinate
            y = self.radius * math.sin(angle)  # Y coordinate
            z = 0.0                            # Z is constant (horizontal plane)
            positions[sensor] = np.array([x, y, z])
            
        print(f"Sensor positions (radius {self.radius} mm, ZX plane):")
        for sensor, pos in positions.items():
            print(f"  Sensor {sensor}: ({pos[0]:.2f}, {pos[1]:.2f}, {pos[2]:.2f}) mm")
            
        return positions

    def load_multi_sensor_data(self, step):
        """Load ray data for all sensors at a specific time step"""
        sensor_data = {}
        
        for sensor in ['A', 'B', 'C']:
            for model in ['1', '2']:
                sensor_name = f"{sensor}{model}"
                filename = f"ray_data_{sensor_name}_step_{step:03d}.txt"
                
                if os.path.exists(filename):
                    try:
                        data = self.load_single_ray_data(filename)
                        if data:
                            sensor_data[sensor_name] = data
                            print(f"Loaded {sensor_name}: {len(data['rays'])} rays")
                    except Exception as e:
                        print(f"Warning: Could not load {filename}: {e}")
                else:
                    print(f"Warning: {filename} not found")
        
        return sensor_data

    def load_single_ray_data(self, filename):
        """Load ray tracing data from a single sensor file"""
        if not os.path.exists(filename):
            return None
            
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
        
        if len(pos_vertices) == 0:
            return None
            
        return {
            'pos_vertices': np.array(pos_vertices),
            'pos_faces': np.array(pos_faces),
            'neg_vertices': np.array(neg_vertices), 
            'neg_faces': np.array(neg_faces),
            'rays': rays
        }

    def apply_sensor_offset(self, vertices, sensor_group):
        """Apply triangular positioning offset to sensor vertices"""
        if sensor_group not in self.sensor_positions:
            return vertices
            
        offset = self.sensor_positions[sensor_group]
        return vertices + offset

    def create_sensor_mesh(self, data, sensor_name, sensor_group):
        """Create trimesh objects for a single sensor with positioning"""
        meshes = []
        ray_paths = []
        
        # Apply sensor positioning offset
        pos_verts_offset = self.apply_sensor_offset(data['pos_vertices'], sensor_group)
        neg_verts_offset = self.apply_sensor_offset(data['neg_vertices'], sensor_group)
        
        # Color scheme for individual sensors
        colors = {
            'A1': [255, 0, 255],   # Magenta
            'A2': [0, 255, 255],   # Cyan
            'B1': [255, 255, 0],   # Yellow
            'B2': [0, 255, 0],     # Green
            'C1': [0, 0, 255],     # Blue
            'C2': [255, 0, 0],     # Red
        }
        
        sensor_color = colors.get(sensor_name, [128, 128, 128])  # Default gray
        
        # Create positive mesh
        if len(data['pos_faces']) > 0:
            pos_mesh = trimesh.Trimesh(vertices=pos_verts_offset, faces=data['pos_faces'])
            pos_colors = np.full((len(data['pos_faces']), 4), 
                               [*sensor_color, 180], dtype=np.uint8)
            pos_mesh.visual.face_colors = pos_colors
            meshes.append(('positive', pos_mesh))
        
        # Create negative mesh (stationary, only add once per sensor group) - GREY
        if len(data['neg_faces']) > 0 and sensor_name.endswith('1'):  # Only add for first model
            neg_mesh = trimesh.Trimesh(vertices=neg_verts_offset, faces=data['neg_faces'])
            neg_colors = np.full((len(data['neg_faces']), 4),
                               [128, 128, 128, 150], dtype=np.uint8)  # Grey
            neg_mesh.visual.face_colors = neg_colors
            meshes.append(('negative', neg_mesh))
        
        # Create rays with offset - use same color as sensor
        hits = 0
        misses = 0
        ray_step = max(1, len(data['rays']) // 200)  # Limit rays for performance
        
        for i, ray in enumerate(data['rays'][::ray_step]):
            is_hit = ray[6] > 0
            
            if is_hit:
                hits += 1
            else:
                misses += 1
                
            # Apply offset to ray positions
            origin = np.array(ray[0:3]) + self.sensor_positions[sensor_group]
            endpoint = np.array(ray[3:6]) + self.sensor_positions[sensor_group]
            
            if is_hit:
                path = trimesh.load_path([origin, endpoint])
                # Use sensor color for hit rays
                path.colors = [[*sensor_color, 200]]
                ray_paths.append(path)
            else:
                # Short orange line for misses
                direction = endpoint - origin
                norm = np.linalg.norm(direction)
                if norm > 0:
                    direction = direction / norm * 2.0
                endpoint_short = origin + direction
                path = trimesh.load_path([origin, endpoint_short])
                path.colors = [[255, 165, 0, 150]]  # Orange for misses
                ray_paths.append(path)
        
        print(f"  {sensor_name}: {hits} hits, {misses} misses, {len(ray_paths)} rays displayed")
        
        return meshes, ray_paths

    def create_scene_for_step(self, step):
        """Create complete scene for a specific time step"""
        print(f"\nCreating scene for step {step}...")
        
        # Load data for all sensors
        sensor_data = self.load_multi_sensor_data(step)
        
        if not sensor_data:
            print(f"No data found for step {step}")
            return None
        
        # Create scene
        scene = trimesh.Scene()
        
        # Process each sensor
        total_meshes = 0
        total_rays = 0
        
        for sensor_name, data in sensor_data.items():
            sensor_group = sensor_name[0]  # 'A', 'B', or 'C'
            
            meshes, ray_paths = self.create_sensor_mesh(data, sensor_name, sensor_group)
            
            # Add meshes to scene
            for mesh_type, mesh in meshes:
                node_name = f"{sensor_name}_{mesh_type}"
                scene.add_geometry(mesh, node_name=node_name)
                total_meshes += 1
            
            # Add rays to scene  
            for i, path in enumerate(ray_paths):
                scene.add_geometry(path, node_name=f"{sensor_name}_ray_{i}")
                total_rays += 1
        
        print(f"Scene created: {total_meshes} meshes, {total_rays} ray paths")
        
        # Add coordinate axes at origin
        axes = trimesh.creation.axis(origin_size=2, axis_radius=0.5, axis_length=15)
        scene.add_geometry(axes, node_name='axes')
        
        return scene

    def find_available_steps(self):
        """Find all available time steps in the current directory"""
        pattern = "ray_data_*_step_*.txt"
        files = glob.glob(pattern)
        
        steps = set()
        for file in files:
            try:
                # Extract step number from filename like "ray_data_A1_step_000.txt"
                parts = file.split('_')
                for i, part in enumerate(parts):
                    if part == 'step' and i + 1 < len(parts):
                        step_with_ext = parts[i + 1]  # "000.txt"
                        step_num = int(step_with_ext.split('.')[0])  # "000" -> 0
                        steps.add(step_num)
                        break
            except (IndexError, ValueError) as e:
                continue
        
        return sorted(list(steps))

    def animate_steps_smooth(self, step_duration=5.0):
        """Show smooth animation through all available time steps in single window"""
        available_steps = self.find_available_steps()
        
        if not available_steps:
            print("No ray data files found!")
            print("Expected files: ray_data_A1_step_000.txt, ray_data_A2_step_000.txt, etc.")
            return
        
        print(f"Found {len(available_steps)} time steps: {available_steps[:10]}{'...' if len(available_steps) > 10 else ''}")
        print(f"Loading all scenes for smooth animation...")
        
        # Pre-load all scenes
        scenes = {}
        for i, step in enumerate(available_steps):
            print(f"Loading step {step} ({i+1}/{len(available_steps)})...", end='\r')
            scene = self.create_scene_for_step(step)
            if scene:
                scenes[step] = scene
        
        print(f"\nLoaded {len(scenes)} scenes successfully")
        
        if not scenes:
            print("No valid scenes could be created!")
            return
        
        print(f"\nStarting smooth animation:")
        print(f"- Close window between steps to continue")
        print(f"- Press ESC or close window to exit")
        
        # Show animation in sequence
        step_keys = list(scenes.keys())
        
        for step in step_keys:
            print(f"\n--- Step {step} ---")
            print("Close the current window to continue to next step...")
            scenes[step].show()
                
        print("Animation complete!")

    def animate_steps(self):
        """Original step-by-step animation (opens new window for each step)"""
        available_steps = self.find_available_steps()
        
        if not available_steps:
            print("No ray data files found!")
            print("Expected files: ray_data_A1_step_000.txt, ray_data_A2_step_000.txt, etc.")
            return
        
        print(f"Found {len(available_steps)} time steps: {available_steps}")
        
        for step in available_steps:
            print(f"\n{'='*50}")
            print(f"STEP {step}")
            print(f"{'='*50}")
            
            scene = self.create_scene_for_step(step)
            if scene:
                print(f"\nDisplaying step {step}")
                print("Controls:")
                print("- Mouse drag: Rotate view") 
                print("- Scroll: Zoom")
                print("- Close window to continue to next step")
                
                scene.show()
            
            # Small delay between steps
            time.sleep(0.5)

    def show_single_step(self, step=0):
        """Show a single time step"""
        scene = self.create_scene_for_step(step)
        if scene:
            print(f"\nDisplaying step {step}")
            scene.show()
        else:
            print(f"Could not create scene for step {step}")

def main():
    print("Multi-Sensor Capacitor Animated Viewer")
    print("-" * 45)
    
    # Initialize viewer with triangular sensor positioning on ZX plane
    viewer = MultiSensorViewer(radius=26.45)  # 26.45mm radius
    
    # Check command line arguments
    if len(sys.argv) > 1:
        if sys.argv[1] == 'animate':
            viewer.animate_steps_smooth()  # New smooth animation
        elif sys.argv[1] == 'animate_old':
            viewer.animate_steps()  # Old step-by-step animation
        else:
            try:
                step = int(sys.argv[1])
                viewer.show_single_step(step)
            except ValueError:
                print("Usage: python multi_sensor_animated_viewer.py [step_number|animate|animate_old]")
                print("  animate     - Smooth animation (recommended)")
                print("  animate_old - Step-by-step with new windows")
                print("  [number]    - Show specific step")
    else:
        # Default: show step 0
        viewer.show_single_step(0)

if __name__ == "__main__":
    main()