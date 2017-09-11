import cobra.test
import warnings

from cobra.flux_analysis import flux_variability_analysis

oxygen = "o2"

acetate = "ac"
acetaldehyde = "acald"
two_oxoglutarate = "akg"
ethanol = "etoh"
d_fructose = "fru"
fumarate = "fum"
d_glucose = "glc__D"
l_glutamine = "gln__L"
l_glutamate = "glu__L"
d_lactate = "lac__D"
l_malate = "mal__L"
pyruvate = "pyr"
succinate = "succ"
ATP = "atp"
ADP = "adp"
ATPM = "ATPM"
NADH = "NADH"
NADPH = "NADPH"
proton = "h"
threePG = "3PG"
PEP="PEP"
PYR="PYR"
OAA="OAA"
G6P="G6P"
F6P="F6P"
R5P="R5P"
E4P="E4P"
G3P="G3P"
ACCOA="ACCOA"
AKG="AKG"
SUCCOA="SUCCOA"
malic_enzyme_NAD = "ME1"
malic_enzyme_NADPH = "ME2"
fumarate_reductaase = "FRD7"
succinate_dehydrogenase = "SUCDi"
malate_dehydrogenase = "MDH"
NAD_transdehydrogenase = "NADTRHD"
pep_carboxykinase = "PPCK"
pyruvate_kinase = "PYK"
PFK = "PFK"
PFK_A_gene = "b1723"
PFK_B_gene = "b3916"
G6PDH2r_gene = "b1852"
ENO_gene = "b2779"

cytosolic_metabolite_suffix = "_c"

glucose_exchange_realistic_lower_bound = -18.5
organic_exchange_realistic_lower_bound = -20.0
unlimited_reaction_bound = 1000.0
limited_reaction_bound = 0.0

zero_threshold = 1E-8

def get_model_soultion(model):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        solution = model.optimize()
        if (not solution.status == "optimal") or solution.objective_value < zero_threshold:
            solution.objective_value=0.0
        return solution

def setup_oxygen_bounds(experiment, aerobic, exact_oxygen_bound):
    if exact_oxygen_bound is not None:
        experiment.reactions_bounds.append(ExactSubstrateReactionBounds(oxygen, exact_oxygen_bound))
        experiment.name = experiment.name + "exact_oxygen"
    elif aerobic is True:
        experiment.name = experiment.name + "aerobic"
        experiment.reactions_bounds.append(UnlimitedSubstrateReactionBounds(oxygen))
    elif aerobic is False:
        experiment.name = experiment.name + "unaerobic"
        experiment.reactions_bounds.append(UnavailableSubstrateReactionBounds(oxygen))


def setup_substrate_bounds(experiment, substrate, exact_substrate_bounds):
    if substrate is not None:
        if exact_substrate_bounds is None:
            experiment.name = experiment.name + "_realistic_" + substrate
            experiment.reactions_bounds.append(RealisticSubstrateReactionBounds(substrate))
        else:
            experiment.name = experiment.name + "_exact_" + substrate
            experiment.reactions_bounds.append(ExactSubstrateReactionBounds(substrate, exact_substrate_bounds))
        if not substrate == d_glucose:
            experiment.reactions_bounds.append(UnavailableSubstrateReactionBounds(d_glucose))


def setup_drain_target(experiment, drain_target):
    if drain_target is not None:
        experiment.reactions_bounds.append(ExactSubstrateReactionBounds(substrate=d_glucose, bound=-1))
        if drain_target == ATPM:
            experiment.reactions_bounds.append(ReactionBounds(reaction_name=ATPM, lower_bound=0))
            experiment.objective_function_name = ATPM
            experiment.drain_target = None
        else:
            experiment.reactions_bounds.append(ExactReactionBounds(reaction_name=ATPM, bound=0))
            experiment.objective_function_name = drain_target + CofactorAndPrecursorsFBATest.drain_suffix
            experiment.drain_target = drain_target.lower()


class FBAExperiment:
    def __init__(self, name,
                 aerobic=True, exact_oxygen_bound=None,
                 substrate=None, exact_substrate_bound=None,
                 drain_target=None,
                 objective_function_name=""):
        self.name = name
        self.reactions_bounds = list()
        self.objective_function_name = objective_function_name
        self.drain_target = None
        setup_oxygen_bounds(self, aerobic, exact_oxygen_bound)
        setup_substrate_bounds(self, substrate, exact_substrate_bound)
        setup_drain_target(self, drain_target)

    def run(self, model):
        return self.set_model_bounds_and_get_optimal_values(model)

    def set_model_bounds(self, model):
        for reaction_bound in self.reactions_bounds:
            self.set_reaction_bounds(model=model, reaction_bound=reaction_bound)

    def set_model_bounds_and_get_optimal_values(self, model):
        self.set_model_bounds(model)
        if self.objective_function_name:
            model.objective = model.reactions.get_by_id(self.objective_function_name)
            get_model_soultion(model)

    def run_experiment_and_print_result(self, model):
        pass

    def set_reaction_bounds(self, reaction_bound, model):
        reaction = model.reactions.get_by_id(reaction_bound.reaction_name)
        if reaction_bound.lower_bound is not None:
            reaction.lower_bound = reaction_bound.lower_bound
        if reaction_bound.upper_bound is not None:
            reaction.upper_bound = reaction_bound.upper_bound


class SubstrateFBAExperiment(FBAExperiment):
    def __init__(self, substrate, exact_substrate_bound=None, aerobic=None, exact_oxygen_bound=None):
        super(SubstrateFBAExperiment, self).__init__(name="",
                                                     substrate=substrate, exact_substrate_bound=exact_substrate_bound,
                                                     aerobic=aerobic, exact_oxygen_bound=exact_oxygen_bound)

    def run(self, model):
        return super(SubstrateFBAExperiment, self).run(model)

    def run_experiment_and_print_result(self, model):
        opt_solution = self.run(model)
        opt_value = opt_solution.objective_value if opt_solution.status == "optimal" else 0.0
        print(self.name + ": " + str(round(opt_value, 4)) + "/h")


class CofactorAndPrecursorsFBATest(FBAExperiment):
    drain_suffix = "_drain"

    def __init__(self, drain_target, aerobic):
        super(CofactorAndPrecursorsFBATest, self).__init__(
            name=drain_target + CofactorAndPrecursorsFBATest.drain_suffix + "_",
            aerobic=aerobic,
            drain_target=drain_target)

    def run(self, model):
        if self.drain_target:
            drain_reaction = cobra.Reaction(self.drain_target.upper() + CofactorAndPrecursorsFBATest.drain_suffix)
            drain_reaction.name = self.drain_target + " drain reaction"
            if self.drain_target.upper() in [NADH, NADPH]:
                oxidized_target = get_cytosolic_metabolite(model, self.drain_target[:-1])
                reduced_target = get_cytosolic_metabolite(model, self.drain_target)
                released_proton = get_cytosolic_metabolite(model, proton)
                drain_reaction.add_metabolites({reduced_target: -1.0,
                                                oxidized_target: 1.0,
                                                released_proton: 1.0})
            else:
                drain_target_metabolite = get_cytosolic_metabolite(model, self.drain_target)
                drain_reaction.add_metabolites({drain_target_metabolite: -1.0})
            model.add_reaction(drain_reaction)
        opt_solution = super(CofactorAndPrecursorsFBATest, self).run(model)
        return opt_solution

    def run_experiment_and_print_result(self, model):
        opt_solution = self.run(model)
        print(self.name + ": " +
              str(round(opt_solution.objective_value, 4)) + "mol/mol_glucose" +
              "  ATP shadow price: " + str(opt_solution.shadow_prices[ADP + cytosolic_metabolite_suffix]))


class FluxVariabilityExperiment(FBAExperiment):
    def __init__(self, substrate, reactions_to_test, aerobic, exact_substrate_bound):
        super(FluxVariabilityExperiment, self).__init__(name="FVA_", substrate=substrate, aerobic=aerobic,
                                                        exact_substrate_bound=exact_substrate_bound)
        self.reactions_to_test = list(reactions_to_test)

    def run(self, model):
        super(FluxVariabilityExperiment, self).set_model_bounds(model)
        return flux_variability_analysis(model, self.reactions_to_test)

    def run_experiment_and_print_result(self, model):
        fva_result = self.run(model)
        print(fva_result)


class KnockoutExperiment(FBAExperiment):
    def __init__(self, substrate=None, aerobic=None, reactions_to_ko=None, genes_to_ko=None):
        name = ""
        if reactions_to_ko is None:
            reactions_to_ko = []
        else:
            name = name + str(reactions_to_ko)
        if genes_to_ko is None:
            genes_to_ko = []
        else:
            name = name + str(genes_to_ko)
        name = name + "_knockout_"
        super(KnockoutExperiment, self).__init__(name=name,
                                                 substrate=substrate,
                                                 aerobic=aerobic)
        self.reactions_to_ko = list(reactions_to_ko)
        self.genes_to_ko = list(genes_to_ko)

    def run(self, model):
        for reaction in self.reactions_to_ko:
            model.reactions.get_by_id(reaction).knock_out()
        for gene in self.genes_to_ko:
            model.genes.get_by_id(gene).knock_out()
        return super(KnockoutExperiment, self).run(model)

    def run_experiment_and_print_result(self, model):
        opt_solution = self.run(model)
        print(self.name, opt_solution.objective_value)


def get_cytosolic_metabolite(model, metabolite):
    return model.metabolites.get_by_id(metabolite + cytosolic_metabolite_suffix)


class ReactionBounds:
    def __init__(self, reaction_name, lower_bound=None, upper_bound=None):
        self.reaction_name = reaction_name
        self.upper_bound = upper_bound
        self.lower_bound = lower_bound


class ExactReactionBounds(ReactionBounds):
    def __init__(self, reaction_name, bound):
        super(ExactReactionBounds, self).__init__(reaction_name, lower_bound=bound, upper_bound=bound)


class RealisticSubstrateReactionBounds(ReactionBounds):
    def __init__(self, substrate):
        super(RealisticSubstrateReactionBounds, self).__init__(exchange_reaction_name_for_substrate(substrate),
                                                               lower_bound=organic_exchange_realistic_lower_bound if not substrate == d_glucose else glucose_exchange_realistic_lower_bound,
                                                               upper_bound=unlimited_reaction_bound)


class UnlimitedSubstrateReactionBounds(ReactionBounds):
    def __init__(self, substrate):
        super(UnlimitedSubstrateReactionBounds, self).__init__(exchange_reaction_name_for_substrate(substrate),
                                                               lower_bound=-unlimited_reaction_bound,
                                                               upper_bound=unlimited_reaction_bound)


class UnavailableSubstrateReactionBounds(ReactionBounds):
    def __init__(self, substrate):
        super(UnavailableSubstrateReactionBounds, self).__init__(exchange_reaction_name_for_substrate(substrate),
                                                                 lower_bound=limited_reaction_bound,
                                                                 upper_bound=unlimited_reaction_bound)


class ExactSubstrateReactionBounds(ExactReactionBounds):
    def __init__(self, substrate, bound):
        super(ExactSubstrateReactionBounds, self).__init__(exchange_reaction_name_for_substrate(substrate),
                                                           bound=bound)


def exchange_reaction_name_for_substrate(substrate_name):
    return "EX_" + substrate_name + "_e"
