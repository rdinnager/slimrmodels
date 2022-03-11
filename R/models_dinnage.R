#' Run a SLiM model with fixed population size dynamics, in multiple connected
#' subpopulations, and return the results.
#'
#' @param n_pop Number of subpopulations to simulate
#' @param n_gen Number of generations to simulate for
#' @param mut_rate Mutation rate
#' @param selection Selection strength
#' @param genome_size Genome size
#' @param recomb_rate Recombination rate
#' @param pop_abund vector or matrix of fixed population sizes. If a matrix, it
#' should be a `n_pop` by `n_gen` matrix giving fixed population sizes of each
#' subpopulation at each generation. If a vector of length one, it will be
#' treated as a constant across all subpopulations and generations. If a vector
#' of length equal to `n_pop` or `n_gen`, it will be assumed to be constant across
#' the other dimension (unless `n_pop` equals `n_gen`, in which case an error
#' will be thrown). Can also be a function that takes the generation as
#' its first argument and returns a vector of population sizes of length `n_pop`
#' @param sampler vector or matrix or function of the same kind described
#' for `pop_abund`, but the values represent the number of individuals to
#' sample for the data output.
#' @param migration_rates An `n_pop` by `n_pop` matrix with fixed
#' @param ... further arguments passed to any functions passed
#' as arguments to `pop_abund`, `sampler`, or `migration_rates`
#'
#' @return
#' @export
#'
#' @examples
mod_fixed_pop_dyn <- function(n_pop = 3, n_gen = 1000, mut_rate = 1e-6, selection = 1e-12, genome_size = 300000,
                              recomb_rate = 1e-8, pop_abund = c(25, 50, 100), sampler = 1,
                              migration_rates = 0.5, ...) {

  slim_script(

    slim_block(initialize(), {

      #setSeed(12345); don't want to set the seed because we want diff each time
      #initializeSLiMOptions(keepPedigrees=T);
      initializeMutationRate(slimr_template("mut_rate", 1e-6));
      initializeMutationType("m1", 0.5, "n", 0, slimr_template("selection_strength", 0.1));
      initializeGenomicElementType("g1", m1, 1.0);
      initializeGenomicElement(g1, 0, slimr_template("genome_size", !!default_genome_size) - 1);
      initializeRecombinationRate(slimr_template("recomb_rate", 1e-8));
      initializeSex("A");
      defineConstant("abund", slimr_inline(pop_abunds, delay = FALSE));
      defineConstant("sample_these", slimr_inline(sample_these, delay = FALSE));
      defineConstant("samp_sizes", slimr_inline(samp_sizes, delay = FALSE));

    }),
    slim_block(1, {

      init_pop = slimr_inline(init_popsize, delay = TRUE)

      ## set populations to initial size
      sim.addSubpop("p1", asInteger(init_pop[0]));
      sim.addSubpop("p2", asInteger(init_pop[1]));
      sim.addSubpop("p3", asInteger(init_pop[2]));

    }),

    slim_block(1, late(), {
      ## get starting population from a file which we will fill-in later
      sim.readFromPopulationFile(slimr_inline(starting_pop, delay = TRUE));
      ## migration on or off flags for pops 1-3 (using tag)
      p1.tag = 0;
      p2.tag = 0;
      p3.tag = 0;
    }),

    slim_block(1, !!n_gen, late(), {

      ## update generation number
      gen = sim.generation %% 50
      if(gen == 0) {
        gen = 50
      }

      ## set population size to observed levels
      p1.setSubpopulationSize(asInteger(ceil(abund[0, gen - 1] * slimr_template("popsize_scaling", 100))));
      p2.setSubpopulationSize(asInteger(ceil(abund[1, gen - 1] * ..popsize_scaling..)));
      p3.setSubpopulationSize(asInteger(ceil(abund[2, gen - 1] * ..popsize_scaling..)));

      ## increase migration when above abundance threshold
      if(p1.tag == 0 & abund[0, gen - 1] > slimr_template("abund_threshold", 5)) {
        p2.setMigrationRates(p1, slimr_template("migration_rate", 0))
        p3.setMigrationRates(p1, ..migration_rate..)
        p1.tag = 1;
      }
      if(p1.tag == 1 & abund[0, gen - 1] <= ..abund_threshold..) {
        p2.setMigrationRates(p1, 0)
        p3.setMigrationRates(p1, 0)
        p1.tag = 0;
      }

      if(p2.tag == 0 & abund[1, gen - 1] > ..abund_threshold..) {
        p1.setMigrationRates(p2, ..migration_rate..)
        p3.setMigrationRates(p2, ..migration_rate..)
        p2.tag = 1;
      }
      if(p2.tag == 1 & abund[1, gen - 1] <= ..abund_threshold..) {
        p1.setMigrationRates(p2, 0)
        p3.setMigrationRates(p2, 0)
        p2.tag = 0;
      }

      if(p3.tag == 0 & abund[2, gen - 1] > ..abund_threshold..) {
        p1.setMigrationRates(p3, ..migration_rate..)
        p2.setMigrationRates(p3, ..migration_rate..)
        p3.tag = 1;
      }
      if(p3.tag == 1 & abund[2, gen - 1] <= ..abund_threshold..) {
        p1.setMigrationRates(p3, 0)
        p2.setMigrationRates(p3, 0)
        p3.tag = 0;
      }

      ## only run if the generation is in our sample_these list
      if(any(match(sample_these, sim.generation) >= 0)) {
        ## find the sample size that matches the matching "year" for our obs data
        ssizes = drop(samp_sizes[ , which(sample_these == sim.generation)])
        ## sample individuals
        ind_sample = c(sample(p1.individuals, ssizes[0]),
                       sample(p2.individuals, ssizes[1]),
                       sample(p3.individuals, ssizes[2]))
        ## output individuals genomes
        slimr_output(ind_sample.genomes.output(), "pop_sample", do_every = 1);
        slimr_output(ind_sample.genomes.individual.subpopulation, "subpops", do_every = 1)
      }

    }),

    slim_block(1000, late(), {
      sim.simulationFinished()
    })

  ) -> pop_sim_samp

}

# 110 line of code originally
