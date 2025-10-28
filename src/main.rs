use macroquad::{prelude::{scene::camera_pos, *}, rand::gen_range, };

//use std::collections::HashMap;
use rustc_hash::FxHashMap as HashMap; //schnellere Hashmap

struct Particle {
    pos: Vec3,
    old_pos: Vec3,
    acceleration: Vec3,
    radius: f32
}
impl Particle {
    fn new(pos: Vec3, velocity: Vec3, radius: f32) -> Self {
        Self {
            pos: pos,
            old_pos: pos - velocity, // alte Position geändert für eine Startgeschwindigkeit
            acceleration: Vec3::ZERO,
            radius,
        }
    }
    fn apply_foce(&mut self, f: Vec3) {
        self.acceleration += f;
    }

    fn update(&mut self, dt: f32, gravity: bool) {
        // Gravitationsbeschleunigung
        if gravity { self.apply_foce(vec3(0., -100., 0.)); } else { self.apply_foce(vec3(0., 0., 0.)); } // nur Gravitation
        // Beschleunigung berechnen. Wird errechnet aus dem Unterschied der alten und aktuellen Position.
        let velocity = self.pos - self.old_pos;

        // Verlet Formel
        let new_pos = self.pos + velocity + self.acceleration * dt * dt; // dt^2 da acc = m / s^2

        // Positionen aktuallisiern.
        self.old_pos = self.pos;

        self.pos = new_pos;

        self.acceleration = vec3(0., 0., 0.)
    }
    
    fn box_constrains(&mut self, size: f32) {
        let damping = 0.5;
        // unten
        if self.pos.y + self.radius > size {
            self.pos.y = size - self.radius;
            self.old_pos.y = self.pos.y + (self.pos.y - self.old_pos.y) * damping;
        }
        // oben
        if self.pos.y - self.radius < 0.0 {
            self.pos.y = self.radius;
            self.old_pos.y = self.pos.y + (self.pos.y - self.old_pos.y) * damping;
        }
        // rechts
        if self.pos.x + self.radius > size {
            self.pos.x = size - self.radius;
            self.old_pos.x = self.pos.x + (self.pos.x - self.old_pos.x) * damping;
        }
        // links
        if self.pos.x - self.radius < 0.0 {
            self.pos.x = self.radius;
            self.old_pos.x = self.pos.x + (self.pos.x - self.old_pos.x) * damping;
        }
        
        // Z achse
        // vorne
        if self.pos.z + self.radius > size {
            self.pos.z = size - self.radius;
            self.old_pos.z = self.pos.z + (self.pos.z - self.old_pos.z) * damping;
        }
        // hinten
        if self.pos.z - self.radius < 0.0 {
            self.pos.z = self.radius;
            self.old_pos.z = self.pos.z + (self.pos.z - self.old_pos.z) * damping;
        }
    }

    fn resolve_collision(&mut self, other: &mut Particle) { // muss neu
        let delta = self.pos - other.pos; // Wie weit sind die beiden auseinander?

        let dist_squared = delta.length_squared(); // errechnet die Länge des Vektors

        let min_dist = self.radius + other.radius;
        let min_dist_squared = min_dist * min_dist;

        if dist_squared - min_dist_squared < 0.0 && dist_squared > 0.0 {
            let dist = dist_squared.sqrt();
            
            let overlap = min_dist - dist;

            let direction = delta / dist;

            let correction = direction * overlap * 0.5;

            self.pos += correction;

            other.pos -= correction;

        }
    }
    fn draw(&mut self, sun_positon: Vec3, camera_pos: Vec3) {
        let velocity = self.pos - self.old_pos;
        let speed_squared = velocity.length_squared();  // ← Nur x² + y²

        let mut color = speed_to_color(speed_squared);

        // Lichteffekt von der Sonnenentferung / Distanzbeleuchtung
        let falloff = 0.008;

        let delta = self.pos - sun_positon;
        let dist = delta.length();

        let mut intensity: f32 = 1.0 / (1.0 + dist * falloff);

        // Rim Lighting
        let to_light = (sun_positon - self.pos).normalize(); // Wo ist das licht aus der Sicht des Partikels?
        let to_camera = (camera_pos - self.pos).normalize(); // Wo ist die Kamera aus der Sicht des Partikels?

        let dot = to_light.dot(to_camera); // Winkel zwischen Licht und Kamera

        let rim = (1.0 - dot.abs()).powf(2.0);

        let rim_boost = 1.0 + rim * 2.; // Heller am Rand

        intensity *= rim_boost;

        // Intensität mit der Farbe multiplizieren
        color.r *= intensity;
        color.g *= intensity;
        color.b *= intensity;
        // color.a *= intensity;


        draw_sphere(
            self.pos, self.radius, None, color
        );
        // draw_sphere_wires(self.pos, self.radius, None, RED);
    }
    fn calm(&mut self) {
        self.old_pos = self.pos;
    }
    
}
// veraltet aber lustig weil es viel weniger Performant ist
fn resolve_collision(particles: &mut Vec<Particle>) {
    for i in 0..particles.len() {
        for j in i+1..particles.len() {
            let (left, right) = particles.split_at_mut(j);
            left[i].resolve_collision(&mut right[0]);
        }
    }
}

fn fill_grid<'a>(hashie: &'a mut HashMap<(i32, i32, i32), Vec<usize>>, particles: &Vec<Particle>, cell_size: f32) -> &'a mut HashMap<(i32, i32, i32), Vec<usize>> {
    //let mut hashie: HashMap<(i32, i32), Vec<usize>> = HashMap::default();
    hashie.reserve(particles.len() / 4);
    for (index, particle) in particles.iter().enumerate() {
        let cell_x = (particle.pos.x / cell_size) as i32;
        let cell_y = (particle.pos.y / cell_size) as i32;
        let cell_z = (particle.pos.z / cell_size) as i32;

        let cell_key = (cell_x, cell_y, cell_z);

        hashie.entry(cell_key)
            .or_insert_with(|| Vec::with_capacity(8))
            .push(index);
    }
    hashie
}
// muss neu das ist der aufwendige part :/
fn resolve_collision_with_grid(particles: &mut Vec<Particle>, grid: &HashMap<(i32, i32, i32), Vec<usize>>) {
    for (cell_key, cell_indices) in grid.iter() { // Über jede Celle wird iteriert.
        for i in 0..cell_indices.len() { // Über alle Particel in einer Cell werden iteriert und verglichen.
            for j in i+1..cell_indices.len() {
                let idx_a = cell_indices[i];  // Konvertiert denn cellindex in einen echten index.
                let idx_b = cell_indices[j];  
                // Finde den größeren und kleineren Index
                let (small_idx, large_idx) = if idx_a < idx_b {
                    (idx_a, idx_b)
                } else {
                    (idx_b, idx_a)
                };
                
                let (left, right) = particles.split_at_mut(large_idx);
                left[small_idx].resolve_collision(&mut right[0]);
            }
        }
        let neighbors = [
            // Gleiche Y-Ebene (8 Nachbarn)
            (cell_key.0 + 1, cell_key.1,     cell_key.2),     // rechts
            (cell_key.0 + 1, cell_key.1,     cell_key.2 + 1), // rechts vorne
            (cell_key.0,     cell_key.1,     cell_key.2 + 1), // vorne
            (cell_key.0 - 1, cell_key.1,     cell_key.2 + 1), // links vorne
            (cell_key.0 - 1, cell_key.1,     cell_key.2),     // links
            (cell_key.0 - 1, cell_key.1,     cell_key.2 - 1), // links hinten
            (cell_key.0,     cell_key.1,     cell_key.2 - 1), // hinten
            (cell_key.0 + 1, cell_key.1,     cell_key.2 - 1), // rechts hinten
            
            // Obere Y-Ebene (9 Nachbarn)
            (cell_key.0,     cell_key.1 + 1, cell_key.2),     // oben mitte
            (cell_key.0 + 1, cell_key.1 + 1, cell_key.2),     // oben rechts
            (cell_key.0 + 1, cell_key.1 + 1, cell_key.2 + 1), // oben rechts vorne
            (cell_key.0,     cell_key.1 + 1, cell_key.2 + 1), // oben vorne
            (cell_key.0 - 1, cell_key.1 + 1, cell_key.2 + 1), // oben links vorne
            (cell_key.0 - 1, cell_key.1 + 1, cell_key.2),     // oben links
            (cell_key.0 - 1, cell_key.1 + 1, cell_key.2 - 1), // oben links hinten
            (cell_key.0,     cell_key.1 + 1, cell_key.2 - 1), // oben hinten
            (cell_key.0 + 1, cell_key.1 + 1, cell_key.2 - 1), // oben rechts hinten
            
            // Untere Y-Ebene (9 Nachbarn)
            (cell_key.0,     cell_key.1 - 1, cell_key.2),     // unten mitte
            (cell_key.0 + 1, cell_key.1 - 1, cell_key.2),     // unten rechts
            (cell_key.0 + 1, cell_key.1 - 1, cell_key.2 + 1), // unten rechts vorne
            (cell_key.0,     cell_key.1 - 1, cell_key.2 + 1), // unten vorne
            (cell_key.0 - 1, cell_key.1 - 1, cell_key.2 + 1), // unten links vorne
            (cell_key.0 - 1, cell_key.1 - 1, cell_key.2),     // unten links
            (cell_key.0 - 1, cell_key.1 - 1, cell_key.2 - 1), // unten links hinten
            (cell_key.0,     cell_key.1 - 1, cell_key.2 - 1), // unten hinten
            (cell_key.0 + 1, cell_key.1 - 1, cell_key.2 - 1), // unten rechts hinten
        ];
        for neighbor_key in neighbors.iter() {
            if let Some(neighbor_indicies) = grid.get(neighbor_key) { // wenn die Nachbarzelle existiert
                for &idx_a in cell_indices.iter() {
                    for &idx_b in neighbor_indicies.iter() {
                        if idx_a == idx_b { continue; } // wenn es die gleichen sind: überspringen
                        let (small_idx, large_idx) = if idx_a < idx_b {
                            (idx_a, idx_b)
                        } else {
                            (idx_b, idx_a)
                        };
                        
                        let (left, right) = particles.split_at_mut(large_idx);
                        left[small_idx].resolve_collision(&mut right[0]);
                    }
                }
            }
        }
    }
}
// passt so
fn update_particles(particles: &mut Vec<Particle>, dt: f32, bool_gravity: bool, box_size: f32) {
    for particle in particles {
        particle.update(dt, bool_gravity);
        particle.box_constrains(box_size);
    }
}

fn draw_particles(particles: &mut Vec<Particle>, sun_position: Vec3, camera_pos: Vec3) {
    for particle in particles {
        particle.draw(sun_position, camera_pos)
    }
}
fn speed_to_color(speed: f32) -> Color {
    let normalized = (speed / 1.0).min(1.0);
    let r = normalized.powf(1.5);
    let g = (1.0 - normalized).powf(2.0);
    Color::new(r, g, 1.0 - r, 1.0)
}
fn camera_movement(
    camera_angle_h: &mut f32,      // ← &mut = mutable Referenz
    camera_angle_v: &mut f32, 
    camera_distance: &mut f32, 
    camera: &mut Camera3D           // ← auch &mut!
) {
    // Jetzt mit * dereferenzieren:
    if is_key_down(KeyCode::Left) {
        *camera_angle_h += 0.02;   // keine ahnung weiso aber naja :/
    }
    if is_key_down(KeyCode::Right) {
        *camera_angle_h -= 0.02;
    }
    if is_key_down(KeyCode::Up) {
        *camera_angle_v += 0.02;
        *camera_angle_v = camera_angle_v.clamp(-1.5, 1.5);
    }
    if is_key_down(KeyCode::Down) {
        *camera_angle_v -= 0.02;
        *camera_angle_v = camera_angle_v.clamp(-1.5, 1.5);
    }
    
    if is_key_down(KeyCode::W) {
        *camera_distance -= 1.0;
    }
    if is_key_down(KeyCode::S) {
        *camera_distance += 1.0;
    }
    *camera_distance = camera_distance.clamp(20.0, 300.0);
    
    let target = vec3(50., 50., 50.);
    camera.position = vec3(
        target.x + *camera_distance * camera_angle_v.cos() * camera_angle_h.sin(),
        target.y + *camera_distance * camera_angle_v.sin(),
        target.z + *camera_distance * camera_angle_v.cos() * camera_angle_h.cos(),
    );
    camera.target = target;
}
fn spawn_particle(box_size: f32) -> Particle { // muss ich mich ran setzen
    let radius = gen_range(2., 5.);
    Particle::new(
        vec3(
            gen_range(radius, box_size - radius),
            gen_range(radius, box_size - radius),
            gen_range(radius, box_size - radius),
        ),
        vec3(0., 0., 0.),
        radius
    )
}
#[macroquad::main("Verlet Partikel")]
async fn main() {
    let mut particles: Vec<Particle> = Vec::with_capacity(5000);

    let mut fps_history: Vec<f32> = Vec::new();

    let mut bool_gravity = true;

    let mut grid: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::default();

    let mut sun_position = vec3(60., 100., -100.);

    let mut camera: Camera3D = Camera3D {
        position: vec3(0., 20., 80.),
        target: vec3(50., 50., 50.),
        up: vec3(0., 1., 0.),
        fovy: 45.,
        ..Default::default()
    };
    let mut camera_angle_h: f32 = 0.0;  // Rotation um Y-Achse
    let mut camera_angle_v: f32 = 0.3;  // Höhe
    let mut camera_distance: f32 = 100.0;

    const FIXED_DT: f32 = 1.0 / 60.0;
    // MAIN LOOP
    loop {
        clear_background(BLACK);
        // Key Inputs
        if is_key_pressed(KeyCode::R) {
            particles = vec![];
        }
        if is_key_pressed(KeyCode::G) {
            if bool_gravity {
                bool_gravity = false;
            } else {
                bool_gravity = true;
            }
        }
        // Kamera
        camera_movement(&mut camera_angle_h, &mut camera_angle_v, &mut camera_distance, &mut camera);

        set_camera(&camera);
        let box_size: f32 = 200.;
        draw_cube_wires(
            vec3(box_size/2., box_size/2., box_size/2.),
            vec3(box_size, box_size, box_size),
            WHITE
        );
        // Light
        if is_key_down(KeyCode::E) {
            sun_position.z += 1.; 
        }
        if is_key_down(KeyCode::D) {
            sun_position.z -= 1.; 
        }
        draw_sphere(sun_position, 10., None,WHITE); // Sonne
        // Physik
        if is_key_pressed(KeyCode::Key1) {
            for _i in 0..10 {
                particles.push(spawn_particle(box_size))
            }
        }
        let substeps = if particles.len() > 20 { 2 } else { 4 }; // mehere Substeps für mehr Stabilität. Weniger Substeps für bessere performance bei hoher Partikelanzahl.
        let sub_dt = FIXED_DT / substeps as f32;
        for _ in 0..substeps {
            update_particles(&mut particles, sub_dt, bool_gravity, box_size);
            
            grid.clear();
            let cell_size = 10.0;
            fill_grid(&mut grid, &particles, cell_size);
            
            resolve_collision_with_grid(&mut particles, &grid);
        }
        // resolve_collision(&mut particles);
        
        draw_particles(&mut particles, sun_position, camera.position);
        
        set_default_camera();
        // HUD
        let fps = get_fps() as f32;
        fps_history.push(fps);
        if fps_history.len() > 30 {
            fps_history.remove(0);
        }
        let avg_fps = fps_history.iter().sum::<f32>() / fps_history.len() as f32;
        let fps_color = if avg_fps > 55.0 { GREEN } else if avg_fps > 30.0 { YELLOW } else { RED };
        draw_text(
            &format!("FPS: {:.1} | Partikel: {}", avg_fps, particles.len()),
            10.0, 20.0, 24.0, fps_color
        );
        next_frame().await
    }
}