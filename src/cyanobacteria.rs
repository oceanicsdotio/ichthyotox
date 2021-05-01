pub mod cyanobacteria {

    
    struct Carbon {

    }

    impl Carbon {

        fn excrete() {

        }

        fn respire() {

        }

        fn synthesize() {

        }

        fn fix() {

        }
    }


    struct Toxin {

    }

    struct LagTox {
        mclr_production_rate: f32,
        mclr_excretion_rate: f32,
        radius: Vec<f32>, // meters
        irradiance: Vec<f32>, // watts per sq meter
        biomass: Vec<f32>, // grams
        delta_rho: Vec<f32>,
        pub carbohydrate: Vec<f32>, // grams
        pub protein: Vec<f32>, // grams
        pub microcystin: Vec<f32>, // grams
        colony_base_radius: f32, // meters
        temp_ref: f32, // reference for limit fcn
        temp_opt: f32, // optimal for growth
        temp_lethal: f32,
        excretion_frac: f32, // unitless
        fixation_max: f32, // per hour rate
        fixation_beta: f32, // shape factor
        respiration_basic: f32, // per hour
        respiration_active: f32, // unitless
        density_max: f32, // empirical kg/m3
        density_min: f32, 
        vesicle_density: f32, // kg/m3
        cell_frac: f32, // volume composed of cell material
        carbon_ratio_max: f32, // empirical
        vesicle_frac: f32, // fraction occupied by vesicles
        irrad_opt: f32, // for grwoth, W/M2
        synthesis_max: f32, // per hour 
        temp_fcn_alpha: f32, // unitless shape factor
        temp_fcn_beta: f32, // shape coef
        cell_density_coefficient: f32, // shape coef
        light_extinction_biomass: f32, // attenuation coef
        light_attenuation_water: f32, // attenuation coef
        shading_upscale: f32
    }




    impl LagTox {
        fn new(
            mclr_excretion_rate: f32,
            mclr_production_rate: f32,
            count: usize,
        ) -> LagTox {
            LagTox{
                mclr_production_rate,
                mclr_excretion_rate,
                radius: vec![0.0; count],
                irradiance: vec![0.0; count],
                biomass: vec![0.0; count],
                carbohydrate: vec![0.0; count],
                protein: vec![0.0; count],
                microcystin: vec![0.0; count],
                delta_rho: vec![0.0; count],
                colony_base_radius: 75e-6,
                temp_ref: 25.0,
                temp_opt: 28.0,
                temp_lethal: 35.0,
                excretion_frac: 0.1,
                fixation_max: 11.4,
                fixation_beta: 2e-2,
                respiration_basic: 4e-3,
                respiration_active: 2e-1,
                density_max: 1150.0,
                density_min: 1037.0,
                vesicle_density: 150.0,
                cell_frac: 0.25,
                carbon_ratio_max: 4.0,
                vesicle_frac: 0.08,
                irrad_opt: 250.0,
                synthesis_max: 0.05,
                temp_fcn_alpha: 0.286,
                temp_fcn_beta: 0.05,
                cell_density_coefficient: 0.7,
                light_extinction_biomass: 14.0,
                light_attenuation_water: 0.15,
                shading_upscale: 1.0
            }
        }

        /**
         * Read count and initial position and state from file
         */
        fn read() {

        }

        /**
         * Write state and position back to file
         */
        fn write() {

        }

        /**
         * Placeholder for previous quick sort implmentation
         */
        fn sort(abs_depth: Vec<f32>, order: Vec<usize>) {

        }

        /**
         * Sorted by depth, starting at surface.

         This is a stencil
         */
        fn carbon_fixation(&mut self, mesh_area: f32, count: usize, global_irradiance: f32) -> Vec<f32> {

            let carbon_fixation: Vec<f32>;
            let irrad_ratio: Vec<f32>;
            let fixation_coef: Vec<f32>;
            let proxy_depth: Vec<f32>;
            let mut avg_self_shade: Vec<f32>;

            for ii in 0..count {
                self.biomass[ii] = (self.carbohydrate[ii] + self.protein[ii]) / mesh_area;

                let effective_mass = -self.light_extinction_biomass*self.biomass[ii];

                avg_self_shade[ii] = effective_mass.exp()-1.0)/effective_mass;

                if (ii == 0) {
                    self.irradiance[ii] = global_irradiance;
                } else {
                    self.irradiance[ii] = self.irradiance[ii-1] * (-light_extinction_biomass*self.shading_upscale*self.biomass[ii-1]).exp();
                }
            }

            for ii in 0..count {

                self.irradiance[ii] *= avg_self_shade[ii] * (self.zpt[ii] * self.light_attenuation_water).exp()

                irrad_ratio[ii] = self.irradiance[ii] / self.irrad_opt;

                fixation_coef[ii] = (2.0 + self.fixation_beta)*irrad_ratio[ii]/(irrad_ratio[ii].powi(2) + self.fixation_beta*irrad_ratio[ii] + 1.0);

                carbon_fixation[ii] = self.fixation_max*fixation_coef[ii]*self.protein[ii]*(1.0 - self.vesicle_frac)*(self.carbon_ratio_max - self.carbohydrate[ii]/self.protein[ii])/self.carbon_ratio_max;

            }

            carbon_fixation
        }

        /**
         * Transfer carbon from external system into the
         * Colony agent. 
         */
        fn carbon_synthesis(self) -> f32 {
            self.carbohydrate * self.synthesis_max * self.temp_limit()
        }

        /**
         * Update protein and dissolved pools due to excretion
         */
        fn carbon_excretion(self) -> f32 {
            self.excretion_frac * self.temp_function() * (self.respiration_basic * self.carbohydrate + self.synthesis_max * self.protein)
        }

        /**
         * Update carbohydrate and dissolved pools due to 
         * respiration
         */
        fn carbon_respiration(self) -> f32 {
            self.respiration_basic * self.temp_function() * self.protein + self.respiration_active*self.synthesis_max*self.temp_limit()*self.carbohydrate
        }

        /**
         * Temperature limitations coefficient in (0,1) for synthesis
         */
        fn temp_limit(self, temperature: f32) -> Vec<f32> {
            (temperature / self.temp_opt * (((temperature - self.temp_lethal)/(self.temp_opt - self.temp_lethal)).powf((self.temp_ref - self.temp_opt) / self.temp_opt ))).powf(4.0)
        }

        fn temp_function(self, temperature: f32) -> f32 {
            self.temp_fcn_alpha * (self.temp_fcn_beta * (temperature - self.temp_opt + self.temp_ref)).exp()
        }

        /**
         * Instantaneous rate of toxin production
         */
        fn microcystin_production(self) -> f32 {
            self.protein * self.mclr_production_rate
        }

        /**
         * Temperature depenedent toxin loss to water column
         */
        fn microcystin_excretion(self) -> f32 {
            self.protein * self.mclr_excretion_rate
        }

        /**
         * Calculate movement due to buoyancy
         * Calls: velocity, zinterp, zlocate, sigma
         */
        fn vertical_movement() {

        }

        /**
         * Random movements due ot turbulence and diffusion
         */
        fn random_walk() {

        }

        /**
         * Velocity of particle in water m/hr
         * Result if positive if lighter than water
         * Calls desnity() and viscoty()
         */
        fn stokes_velocity(self) -> f32 {
            60.0*60.0*2.0/9.0*9.81*self.radius.powi(2)*(self.rho-self.density)/self.dynamic_viscosity()
        }

        fn algae_density(self) -> f32 {
            let cell_density: f32 = self.density_min + (self.density_max - self.density_min)*(1.0 - (-self.cell_density_coefficient * self.carbohydrate / self.protein).exp());

            (1.0 - self.cell_frac)*(self.rho + 0.7) + self.cell_frac*((1.0 - self.vesicle_frac)*cell_density + self.vesicle_frac * self.vesicle_density)

        }

        fn dynamic_viscosity(self, simple: bool, temperature: f32) {

            let mut viscosity: f32;
            if (simple) {
                viscosity = 1e-3 * 10.0f32.powf(-1.65 + 262.0/(temperature + 169.0))
            } else {
                let aa = 1.541 + 19.998*0.01*temperature - 9.52*10.0f32.powi(-5) * temperature.powi(2);

                let bb = 7.974 - 7.561*0.01 + 4.724*0.0001 * temperature.powi(2);

                let vsic_pure = 4.2844e-5 + 1.0/(0.157 * (temperature + 64.993).powi(2) - 91.296);

                viscosity = vsic_pure*(1.0 + aa*self.sal + bb*self.sal.powi(2));
            }

            viscosity

        }
    }
}
