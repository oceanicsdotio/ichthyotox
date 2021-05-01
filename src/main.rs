
mod cyanobacteria;
mod behavior;
mod variables;


pub mod MOD_LAG {

    struct LAG_SIM {
        pub simID: String,
        pub nnodes: u8,
        pub nelements: u8,
        pub nlayers: u8,
        pub lines_read: u8,
        pub globalIrradiance: f32,
        pub meshArea: f32,
        pub layerDepth: f32,
        pub layerSigma: f32,
        pub time: f32,
        pub daytime: f32,
        pub clocktime: f32,
        elementSigmaVolume: Vec<f32>,
        elementArea: Vec<f32>,
        pub verticaltox: Vec<f32>,
        pub verticaldiff: Vec<f32>,
        pub verticaltemp: Vec<f32>,
        pub verticalrho: Vec<f32>
    }

    impl LAG_SIM {
        fn init(
            simID: String,
            toxin: f32, 
            nlayers: usize,
            nelements: u8,
            nnodes: u8

        ) -> LAG_SIM {

           

            LAG_SIM {
                simID,
                nnodes,
                nelements,
                nlayers: nlayers as u8,
                lines_read: 0,
                globalIrradiance: 0.0,
                meshArea: 0.0,
                layerDepth: 0.0,
                layerSigma: 0.0,
                time: 0.0,
                daytime: 0.0,
                clocktime: 0.0,
                elementSigmaVolume: Vec::new(),
                elementArea: Vec::new(),
                verticaltox: vec![toxin; nlayers],
                verticaldiff: vec![0.0; nlayers],
                verticaltemp: vec![0.0; nlayers],
                verticalrho: vec![0.0; nlayers]
            }
        }

        /**
         * Read physical forcing conditions from a file:

         - u_vel, v_vel, w_vel, diffusivity, elevation, salinity, tempature, density
         */
        fn read(self) {

        }

        fn diffuse() {

        }


    }

    fn triangle_grid_edge() {

    }

    fn hunt() {

    }

    fn spline() {

    }

}


fn main() {


    struct Experiment {
        temperature: f32,
        slope: f32,
        folder_prefix: String, 
        dt: f32,
        diffusivity: f32
    }



    impl Experiment {

        /**
         * Millero and Poisson
         */
        fn density(
            temperature: f32, 
            salinity: f32
        ) -> f32 {


            let aa = 999.842594 + 6.793952e-2 * temperature - 9.09529e-3 * temperature.powi(2) + 1.001685e-4 * temperature.powi(3) - 1.120083e-6 * temperature.powi(4) + 6.536332e-3 * temperature.powi(5);

            let bb = salinity * (0.824493 - 4.0899e-3*temperature + 7.6438e-5*temperature.powi(2) - 8.2467e-7 * temperature.powi(3) + 5.3875e-9 * temperature.powi(4));

            let cc = salinity.powf(1.5) * -0.00572466 + 0.00010227 * temperature - 1.6546e-6 * temperature.powi(2);

            let dd = 4.8314e-4 * salinity.powi(2);

            aa + bb + cc + dd
        }


        fn new() -> Experiment {
            Experiment {
                temperature: 20.0,
                slope: 694e-5,
                folder_prefix: ".".to_string(), 
                dt: 0.1,
                diffusivity: 3600e-5
            }
        }

        fn temperature(&self, step: usize) -> f32 {
            self.temperature + (step as f32) * self.dt * self.slope
        }

    }

    let experiment = Experiment::new();

    // skip creating mesh


    // diffusivity, temperature, rho_pure, density



    println!("Hello, world!");
}
