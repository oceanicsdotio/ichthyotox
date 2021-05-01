pub mod cyanobacteria {


    struct Coefficients {
        α: f32,
        β: f32
    }
    
    /**
     * Regulate vital rates based on the thermodynamics
     * of the agent. This is generally abstracted to 
     * a saturation curve.
     */
    struct ThermalSystem {
        reference: f32,
        optimal: f32,
        lethal: f32,
        coefficients: Coefficients
    }

    impl ThermalSystem {
        fn new() -> ThermalSystem {
            ThermalSystem {
                reference: 25.0,
                optimal: 28.0,
                lethal: 35.0,
                coefficients: Coefficients {
                    α: 0.286,
                    β: 0.05
                }
            }
        }

        /**
         * Temperature limitations coefficient in (0,1) 
         * for synthesis
         */
        fn limit(&self, temperature: &f32) -> f32 {
            (temperature / self.optimal * (((temperature - self.lethal)/(self.optimal - self.lethal)).powf((self.reference - self.optimal) / self.optimal ))).powf(4.0)
        }

        fn function(&self, temperature: &f32) -> f32 {
            self.coefficients.α * (self.coefficients.β * (temperature - self.optimal + self.reference)).exp()
        }
    }

    struct Buoyancy {
        fraction: f32,
        density: f32
    }

    struct BuoyancySystem {
        radius: f32,
        vesicle: Buoyancy,
        cell: Buoyancy
    }

    impl BuoyancySystem {

        fn new() -> BuoyancySystem {
            BuoyancySystem {
                vesicle: Buoyancy{
                    fraction: 0.08,
                    ρ: 150.0,
                },
                cell: Buoyancy{
                    fraction: 0.25,
                    ρ: 0.0
                }
            }
        }

        fn equilibrium(&self, ρ: f32) -> f32 {
            let eq_density: f32 = ((ρ - (1.0 - self.cell.fraction)*(ρ + 0.7))/self.cell.fraction - self.vesicle_frac(self.vesicle.ρ)/(1.0 - self.vesicle.fraction));

            (1.0 - (eq_density - self.density_min)/(self.density_max - self.density_min)).log()/(-self.cell_density_coefficient)
        }
    }

    



    struct LagTox {
        zpt: f32,
        thermal: ThermalSystem,
        buoyancy: BuoyancySystem,
        mclr_production_rate: f32,
        mclr_excretion_rate: f32,
        radius: f32, // meters
        pub carbohydrate: f32, // grams
        pub protein: f32, // grams
        pub microcystin: f32, // grams
        excretion_frac: f32, // unitless
        fixation_max: f32, // per hour rate
        fixation_beta: f32, // shape factor
        respiration_basic: f32, // per hour
        respiration_active: f32, // unitless
        density_max: f32, // empirical kg/m3
        density_min: f32, 
        carbon_ratio_max: f32, // empirical
        irrad_opt: f32, // for grwoth, W/M2
        synthesis_max: f32, // per hour 
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
            biomass: f32,
            microcystin: f32,
            zpt: f32
        ) -> LagTox {

            let protein = biomass / count;
            let buoyancy = BuoyancySystem::new();

            LagTox{
                zpt,
                thermal: ThermalSystem::new(),
                buoyancy,
                mclr_production_rate,
                mclr_excretion_rate,
                radius: 75e-6,
                carbohydrate: buoyancy.equilibrium()*protein,
                protein,
                microcystin: microcystin / count,
                excretion_frac: 0.1,
                fixation_max: 11.4,
                fixation_beta: 2e-2,
                respiration_basic: 4e-3,
                respiration_active: 2e-1,
                density_max: 1150.0,
                density_min: 1037.0,
                carbon_ratio_max: 4.0,
                irrad_opt: 250.0,
                synthesis_max: 0.05,
                cell_density_coefficient: 0.7,
                light_extinction_biomass: 14.0,
                light_attenuation_water: 0.15,
                shading_upscale: 1.0
            }
        }

        

        fn biomass(&self, area: &f32) -> f32 {
            (self.carbohydrate + self.protein) / area
        }


        /**
         * Stencilself.biomass(area)
         */ 
        fn light_below(
            &self, 
            area: &f32,
            irradiance: &f32
        ) -> f32 {
            irradiance * (-self.light_extinction_biomass * self.shading_upscale * self.biomass(area)).exp()
        }
    
        /**
         * Sorted by depth, starting at surface.

         This is a stencil
         */
        fn carbon_fixation(
            &self, area: &f32, 
            irradiance: &f32
        ) -> f32 {

            let carbon_fixation: f32;
            let fixation_coef: f32;
            let proxy_depth: f32;
            
            let effective_mass: f32 = -self.light_extinction_biomass * (self.biomass(area));

            let irrad_ratio: f32 = irradiance * (effective_mass.exp()-1.0) / effective_mass * (self.zpt * self.light_attenuation_water).exp() / self.irrad_opt;

            let fixation_coef: f32 = (2.0 + self.fixation_beta)*irrad_ratio/(irrad_ratio.powi(2) + self.fixation_beta*irrad_ratio + 1.0);

            self.fixation_max*fixation_coef*self.protein*(1.0 - self.vesicle_frac)*(self.carbon_ratio_max - self.carbohydrate/self.protein)/self.carbon_ratio_max
        }

        /**
         * Transfer carbon from external system into the
         * Colony agent. 
         */
        fn carbon_synthesis(&self, temperature: &f32) -> f32 {
            self.carbohydrate * self.synthesis_max * self.thermal.limit(temperature)
        }

        /**
         * Update protein and dissolved pools due to excretion
         */
        fn carbon_excretion(&self, temperature: &f32) -> f32 {
            self.excretion_frac * self.thermal.function(temperature) * (self.respiration_basic * self.carbohydrate + self.synthesis_max * self.protein)
        }

        /**
         * Update carbohydrate and dissolved pools due to 
         * respiration
         */
        fn carbon_respiration(&self, temperature: &f32) -> f32 {
            self.respiration_basic * self.thermal.function(temperature) * self.protein + self.respiration_active*self.synthesis_max*self.thermal.limit(temperature)*self.carbohydrate
        }

        

        /**
         * Instantaneous rate of toxin production
         */
        fn microcystin_production(&self) -> f32 {
            self.protein * self.mclr_production_rate
        }

        /**
         * Temperature depenedent toxin loss to water column
         */
        fn microcystin_excretion(&self) -> f32 {
            self.protein * self.mclr_excretion_rate
        }

        /**
         * Calculate movement due to buoyancy
         * Calls: velocity, zinterp, zlocate, sigma
         */
        fn vertical_movement(&self) -> f32 {
            0.0
        }

        /**
         * Random chnange in depth due ot turbulence and diffusion. Ross and Sharples 2004
         */
        fn random_walk(&self, dt: &f32, random: &f32) -> f32 {
            let kzp: f32 = 60.0*60.0*0.00001;
            let dkzp: f32 = 0.0;
            let variance: f32 = 1.0;
            
            dkzp * dt + random * ((2.0*kzp + dkzp.powi(2))*dt /variance).sqrt()

        }

        /**
         * Velocity of particle in water m/hr
         * Result if positive if lighter than water
         * Calls desnity() and viscoty()
         */
        fn stokes_velocity(&self, simple: bool, temperature: &f32, salinity: &f32, density: &f32) -> f32 {
            60.0*60.0*2.0/9.0*9.81*self.radius.powi(2)*(density-self.algae_density(density))/self.dynamic_viscosity(simple, temperature, salinity)
        }

        /**
         * Density of cellular material
         */
        fn cell_density(&self) -> f32 {
            self.density_min + (self.density_max - self.density_min)*(1.0 - (-self.cell_density_coefficient * self.carbohydrate / self.protein).exp())
        }

        /**
         * Density of all bloom material
         */ 
        fn algae_density(&self, density: &f32) -> f32 {
            (1.0 - self.cell_frac)*(density + 0.7) + self.cell_frac*((1.0 - self.vesicle_frac)*self.cell_density() + self.vesicle_frac * self.vesicle_density)
        }

        fn dynamic_viscosity(&self, simple: bool, temperature: &f32, salinity: &f32) -> f32 {

            let mut viscosity: f32;
            if (simple) {
                viscosity = 1e-3 * 10.0f32.powf(-1.65 + 262.0/(temperature + 169.0))
            } else {
                let aa = 1.541 + 19.998*0.01*temperature - 9.52*10.0f32.powi(-5) * temperature.powi(2);

                let bb = 7.974 - 7.561*0.01 + 4.724*0.0001 * temperature.powi(2);

                let vsic_pure = 4.2844e-5 + 1.0/(0.157 * (temperature + 64.993).powi(2) - 91.296);

                viscosity = vsic_pure*(1.0 + aa*salinity + bb*salinity.powi(2));
            }

            viscosity

        }
    }
}
