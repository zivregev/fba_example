import fba_example as fba
import cobra.test
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion)

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LightSource


def run_experiments_and_print_result(model, experiments):
    for experiment in experiments:
        with model as model:
            experiment.run_experiment_and_print_result(model)


def run_example_1():
    print("Maximal growth rate on organic substrates:")
    model = cobra.test.create_test_model("textbook")
    substrates = [fba.acetate, fba.acetaldehyde, fba.two_oxoglutarate, fba.ethanol, fba.d_fructose, fba.fumarate,
                  fba.d_glucose, fba.l_glutamine,
                  fba.l_glutamate, fba.d_lactate, fba.l_malate, fba.pyruvate, fba.succinate]
    experiments = []
    for substrate in substrates:
        for aerobic in [True, False]:
            experiments.append(fba.SubstrateFBAExperiment(substrate=substrate, aerobic=aerobic))
    run_experiments_and_print_result(model, experiments)


def run_example_2():
    print("Maximum yield of cofactors from glucose")
    model = cobra.test.create_test_model("textbook")
    experiments = []
    for cofactor in [fba.ATPM, fba.NADH, fba.NADPH]:
        for aerobic in [True, False]:
            experiments.append(fba.CofactorAndPrecursorsFBATest(cofactor, aerobic))
    run_experiments_and_print_result(model, experiments)


def run_example_2_1():
    print("Maximum yield of biosynthetic precursors from glucose")
    model = cobra.test.create_test_model("textbook")
    experiments=[]
    for precursor in [fba.threePG, fba.PEP, fba.PYR, fba.OAA, fba.G6P, fba.F6P, fba.R5P, fba.E4P, fba.G3P, fba.ACCOA, fba.AKG, fba.SUCCOA]:
        experiments.append(fba.CofactorAndPrecursorsFBATest(drain_target=precursor,aerobic=True))
    run_experiments_and_print_result(model,experiments)


def run_example_3_1():
    print("Alternate optimal solutions")
    model = cobra.test.create_test_model("textbook")
    model.reactions.get_by_id("EX_glc__D_e").lower_bound = 0
    model.reactions.get_by_id("EX_succ_e").lower_bound = -20
    model.reactions.get_by_id("EX_succ_e").upper_bound = -20
    solution = model.optimize()
    model.reactions.get_by_id("Biomass_Ecoli_core").lower_bound = solution.f
    model.reactions.get_by_id("Biomass_Ecoli_core").upper_bound = solution.f
    model.objective = model.reactions.get_by_id("ME1")
    min_solution = model.optimize(objective_sense="minimize")
    max_solution = model.optimize(objective_sense="maximize")
    print(min_solution.objective_value)
    print(max_solution.objective_value)


def run_example_3():
    print("Alternate optimal solutions")
    model = cobra.test.create_test_model("textbook")
    reactions_to_test = [model.reactions.get_by_id(reaction) for reaction in {fba.malic_enzyme_NAD,
                                                                              fba.malic_enzyme_NADPH,
                                                                              fba.fumarate_reductaase,
                                                                              fba.succinate_dehydrogenase,
                                                                              fba.malate_dehydrogenase,
                                                                              fba.NAD_transdehydrogenase,
                                                                              fba.pep_carboxykinase,
                                                                              fba.pyruvate_kinase}]
    experiment = fba.FluxVariabilityExperiment(substrate=fba.succinate,
                                               reactions_to_test=reactions_to_test,
                                               aerobic=True,
                                               exact_substrate_bound=-20.0)
    run_experiments_and_print_result(model=model, experiments=[experiment])


def run_example_4_1():
    print("Robustness analysis, variable glucose")
    model = cobra.test.create_test_model("textbook")
    growth_rates = []
    for i in range(0, 21):
        with model as model:
            growth_rates.append(fba.SubstrateFBAExperiment(substrate=fba.d_glucose,
                                                           exact_substrate_bound=-i,
                                                           exact_oxygen_bound=-17)
                                .run(model).objective_value)
    plt.plot(growth_rates)
    plt.show()


def run_example_4_2():
    print("Robustness analysis, variable oxygen")
    model = cobra.test.create_test_model("textbook")
    growth_rates = []
    for i in range(0, 26):
        with model as model:
            growth_rates.append(fba.SubstrateFBAExperiment(substrate=fba.d_glucose,
                                                           exact_substrate_bound=-10,
                                                           exact_oxygen_bound=-i)
                                .run(model).objective_value)
    plt.plot(growth_rates)
    plt.show()


def run_example_5():
    print("Phenotypic phase analysis")
    model = cobra.test.create_test_model("textbook")
    glucose_upper_limit = 21
    oxygen_upper_limit = 21
    glucoses = np.arange(glucose_upper_limit)
    oxygens = np.arange(oxygen_upper_limit)
    glucoses, oxygens = np.meshgrid(glucoses, oxygens)
    growth_rates = calc_growth_rate(glucoses, oxygens, model)

    # TODO pretty graphics
    H = np.histogram(growth_rates, bins=5)

    plt.clf()

    light = LightSource(45, 60)

    white = np.ones((growth_rates.shape[0], growth_rates.shape[1], 3))
    green = white * np.array([0, 1, 0])

    illuminated_surface = light.shade_rgb(green, growth_rates)
    ax = plt.axes(projection='3d')
    ax.plot_surface(glucoses, oxygens, growth_rates, cstride=1, rstride=1, facecolors=illuminated_surface)
    plt.show()
    plt.clf()

    plt.imshow(growth_rates, cmap='hot', interpolation='nearest')
    plt.show()
    plt.clf()


@np.vectorize
def calc_growth_rate(glucose_bound, oxygen_bound, model):
    with model as model:
        fba_solution = fba.SubstrateFBAExperiment(substrate=fba.d_glucose, exact_substrate_bound=-glucose_bound,
                                                  exact_oxygen_bound=-oxygen_bound).run(model)
        return fba_solution.objective_value if fba_solution.status == "optimal" else 0.0


def run_example_6_1():
    model = cobra.test.create_test_model("textbook")
    print("Textbook: ", model.optimize().objective_value)
    run_experiments_and_print_result(model=model,
                                     experiments=[fba.KnockoutExperiment(reactions_to_ko=[fba.PFK], aerobic=True),
                                                  fba.KnockoutExperiment(genes_to_ko=[fba.PFK_A_gene], aerobic=True),
                                                  fba.KnockoutExperiment(genes_to_ko=[fba.PFK_B_gene], aerobic=True),
                                                  fba.KnockoutExperiment(genes_to_ko=[fba.PFK_A_gene, fba.PFK_B_gene], aerobic=True),
                                                  fba.SubstrateFBAExperiment(substrate=fba.d_glucose, aerobic=True),
                                                  fba.KnockoutExperiment(substrate=fba.d_glucose, genes_to_ko=[fba.G6PDH2r_gene], aerobic=True),
                                                  fba.KnockoutExperiment(substrate=fba.d_glucose, genes_to_ko=[fba.ENO_gene], aerobic=True)])


def run_example_6_2():
    model = cobra.test.create_test_model("textbook")
    single_deletion_results = single_gene_deletion(cobra_model=model,gene_list=[fba.ENO_gene,fba.PFK_B_gene])
    print(single_deletion_results)


def run_example_6_3():
    model = cobra.test.create_test_model("textbook")
    double_deletion_results = double_gene_deletion(cobra_model=model, return_frame=True, number_of_processes=1)
    plt.imshow(double_deletion_results)
    plt.colorbar()
    plt.show()


def run_example_7():
    model = cobra.test.create_test_model("textbook")
    biomass_reaction=model.reactions.get_by_id("Biomass_Ecoli_core")
    for metabolite in biomass_reaction.metabolites.keys():
        print(metabolite, ":", biomass_reaction.metabolites.get(metabolite))
    print(biomass_reaction.reaction)
