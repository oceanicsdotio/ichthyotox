




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



    println!("Hello, world!");
}
