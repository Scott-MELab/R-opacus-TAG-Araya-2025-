"""# LIB"""

import cobra
from cobra.io import read_sbml_model, write_sbml_model, save_json_model
from cobra import Model, Reaction, Metabolite
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

model = Model('Core_DSM43205')

"""# Metabolites"""

model.add_metabolites([
    Metabolite('Acetate_c', name='Acetate_c', formula='', compartment='c', charge=0),
    Metabolite('ATP', name='ATP', formula='', compartment='c', charge=0),
    Metabolite('Acetyl_CoA', name='Acetyl_CoA', formula='', compartment='c', charge=0),
    Metabolite('ADP', name='ADP', formula='', compartment='c', charge=0),
    Metabolite('CoA', name='CoA', formula='', compartment='c', charge=0),
    Metabolite('NADH', name='NADH', formula='', compartment='c', charge=0),
    Metabolite('NAD', name='NAD', formula='', compartment='c', charge=0),
    Metabolite('Glycerol_3P', name='Glycerol_3P', formula='', compartment='c', charge=0),
    Metabolite('Proton_c', name='Proton_c', formula='', compartment='c', charge=0),
    Metabolite('Phosphoenolpyruvate', name='Phosphoenolpyruvate', formula='', compartment='c', charge=0),
    Metabolite('Oxaloacetate', name='Oxaloacetate', formula='', compartment='c', charge=0),
    Metabolite('Water_c', name='Water_c', formula='', compartment='c', charge=0),
    Metabolite('Citrate', name='Citrate', formula='', compartment='c', charge=0),
    Metabolite('Cis_Aconitate', name='Cis_Aconitate', formula='', compartment='c', charge=0),
    Metabolite('Isocitrate', name='Isocitrate', formula='', compartment='c', charge=0),
    Metabolite('NADPH', name='NADPH', formula='', compartment='c', charge=0),
    Metabolite('NADP', name='NADP', formula='', compartment='c', charge=0),
    Metabolite('Two_Ketoglutarate', name='Two_Ketoglutarate', formula='', compartment='c', charge=0),
    Metabolite('Phosphate', name='Phosphate', formula='', compartment='c', charge=0),
    Metabolite('Succinyl_CoA', name='Succinyl_CoA', formula='', compartment='c', charge=0),
    Metabolite('Succinate', name='Succinate', formula='', compartment='c', charge=0),
    Metabolite('Ubiquinone8', name='Ubiquinone8', formula='', compartment='m', charge=0),
    Metabolite('Ubiquinol8', name='Ubiquinol8', formula='', compartment='m', charge=0),
    Metabolite('Fumarate', name='Fumarate', formula='', compartment='c', charge=0),
    Metabolite('Malate', name='Malate', formula='', compartment='c', charge=0),
    Metabolite('Glyoxylate', name='Glyoxylate', formula='', compartment='c', charge=0),
    Metabolite('CO2_c', name='CO2_c', formula='', compartment='c', charge=0),
    Metabolite('Propionyl_CoA', name='Propionyl_CoA', formula='', compartment='c', charge=0),
    Metabolite('Proton_e', name='Proton_e', formula='', compartment='e', charge=0),
    Metabolite('Oxygen_c', name='Oxygen_c', formula='', compartment='c', charge=0),
    Metabolite('Diphosphate', name='Diphosphate', formula='', compartment='c', charge=0),
    Metabolite('Oxidized_ferredoxins', name='Oxidized_ferredoxins', formula='', compartment='c', charge=0),
    Metabolite('Reduced_ferredoxins', name='Reduced_ferredoxins', formula='', compartment='c', charge=0),
    Metabolite('Pyruvate', name='Pyruvate', formula='', compartment='c', charge=0),
    Metabolite('AMP', name='AMP', formula='', compartment='c', charge=0),
    Metabolite('Malonyl_CoA', name='Malonyl_CoA', formula='', compartment='c', charge=0),
    Metabolite('Tetradecanoate', name='Tetradecanoate', formula='', compartment='c', charge=0),
    Metabolite('Palmitate', name='Palmitate', formula='', compartment='c', charge=0),
    Metabolite('Stereate', name='Stereate', formula='', compartment='c', charge=0),
    Metabolite('Pentadecanoate', name='Pentadecanoate', formula='', compartment='c', charge=0),
    Metabolite('Heptadecanoate', name='Heptadecanoate', formula='', compartment='c', charge=0),
    Metabolite('Margaroleate', name='Margaroleate', formula='', compartment='c', charge=0),
    Metabolite('Palmitoleate', name='Palmitoleate', formula='', compartment='c', charge=0),
    Metabolite('Cis_Trans_Vaccenate', name='Cis_Trans_Vaccenate', formula='', compartment='c', charge=0),
    Metabolite('Acetate_e', name='Acetate_e', formula='', compartment='e', charge=0),
    Metabolite('Water_e', name='Water_e', formula='', compartment='e', charge=0),
    Metabolite('Oxygen_e', name='Oxygen_e', formula='', compartment='e', charge=0),
    Metabolite('CO2_e', name='CO2_e', formula='', compartment='e', charge=0),
    Metabolite('Glucose_e', name='Glucose_e', formula='', compartment='e', charge=0),
    Metabolite('D_Glucose_6_phosphate', name='D_Glucose_6_phosphate', formula='', compartment='c', charge=0),
    Metabolite('D_Fructose_6_phosphate', name='D_Fructose_6_phosphate', formula='', compartment='c', charge=0),
    Metabolite('D_Fructose_1_6_bisphosphate', name='D_Fructose_1_6_bisphosphate', formula='', compartment='c', charge=0),
    Metabolite('Glyceraldehyde_3_phosphate', name='Glyceraldehyde_3_phosphate', formula='', compartment='c', charge=0),
    Metabolite('Dihydroxyacetone_phosphate', name='Dihydroxyacetone_phosphate', formula='', compartment='c', charge=0),
    Metabolite('Phospho3_D_glyceroyl_phosphate', name='3_Phospho_D_glyceroyl_phosphate', formula='', compartment='c', charge=0),
    Metabolite('Phospho3_D_glycerate', name='3_Phospho_D_glycerate', formula='', compartment='c', charge=0),
    Metabolite('D_Glycerate_2_phosphate', name='D_Glycerate_2_phosphate', formula='', compartment='c', charge=0),
    Metabolite('Pentadecenoic_acid', name = 'Pentadecenoic acid', formula = '', compartment = 'c', charge = 0),
    Metabolite('D_Ribulose_1_5_bisphosphate', name = 'D_Ribulose_1,5_bisphosphate', formula = '', compartment = 'c', charge = 0),
    Metabolite('D_Erythrose_4_phosphate', name = 'D_Erythrose_4_phosphate', formula = '', compartment = 'c', charge = 0),
    Metabolite('Sedoheptulose_1_7_bisphosphate', name = 'Sedoheptulose_1_7_bisphosphate', formula = '', compartment = 'c', charge = 0),
    Metabolite('Sedoheptulose_7_phosphate', name = 'Sedoheptulose_7_phosphate', formula = '', compartment = 'c', charge = 0),
    Metabolite('Alpha_D_Ribose_5_phosphate', name = 'Alpha_D_Ribose_5_phosphate', formula = '', compartment = 'c', charge = 0),
    Metabolite('D_Ribulose_5_phosphate', name = 'D_Ribulose_5_phosphate', formula = '', compartment = 'c', charge = 0),
    Metabolite('D_Xylulose_5_phosphate', name = 'D_Xylulose_5_phosphate', formula = '', compartment = 'c', charge = 0),
    Metabolite('Acetyl_phosphate', name = 'Acetyl phosphate', formula = '', compartment = 'c', charge = 0),
    Metabolite('H2_c', name = 'H2_c', formula = '', compartment = 'c', charge = 0),
    Metabolite('H2_e', name = 'H2_e', formula = '', compartment = 'e', charge = 0),
    Metabolite('Propanol_e', name = 'Propanol_e', formula = '', compartment = 'e', charge = 0),
    Metabolite('Propanol_c', name = 'Propanol_c', formula = '', compartment = 'c', charge = 0),
    Metabolite('Butanoate_c', name = 'Butanoate_c', formula = '', compartment = 'c', charge = 0),
    Metabolite('Butanoate_e', name = 'Butanoate_e', formula = '', compartment = 'e', charge = 0),
    Metabolite('Butanoyl_CoA', name = 'Butanoyl_CoA', formula = '', compartment = 'c', charge = 0),
    Metabolite('Valerate_c', name = 'Valerate_c', formula = '', compartment = 'c', charge = 0),
    Metabolite('Valerate_e', name = 'Valerate_e', formula = '', compartment = 'e', charge = 0),
    Metabolite('Valerly_CoA', name = 'Valerly_CoA', formula = '', compartment = 'c', charge = 0),
    Metabolite('Methanol_e', name = 'Methanol_e', formula = '', compartment = 'e', charge = 0),
    Metabolite('Methanol_c', name = 'Methanol_c', formula = '', compartment = 'c', charge = 0),
    Metabolite('Formaldehyde', name = 'Formaldehyde', formula = '', compartment = 'c', charge = 0),
    Metabolite('D_glucono_1_5_lactone_6_phosphate', name = '6-phospho D-glucono-1,5-lactone	', formula = '', compartment = 'c', charge = 0),
    Metabolite('D_gluconate_6_phosphate', name = 'D_gluconate_6_phosphate', formula = '', compartment = 'c', charge = 0),
    Metabolite('KDPG', name = '2-dehydro-3-deoxy-D-gluconate 6-phosphate', formula = '', compartment = 'c', charge = 0),
    Metabolite('GTP', name = 'GTP', formula = '', compartment = 'c', charge = 0),
    Metabolite('GDP', name = 'GDP', formula = '', compartment = 'c', charge = 0),
    Metabolite('Hydrogen_carbonate', name = 'Hydrogen carbonate', formula = '', compartment = 'c', charge = 0),
    Metabolite('Glucose_1_phosphate', name = 'glucose 1-phosphate', formula = '', compartment = 'c', charge = 0),
    Metabolite('ADPglucose', name = 'ADPglucose', formula = '', compartment = 'c', charge = 0),
    Metabolite('Glycogen_e', name = 'Glycogen_e', formula = '', compartment = 'e', charge = 0),
    Metabolite('Glycogen_c', name = 'Glycogen_c', formula = '', compartment = 'c', charge = 0),

])

for metabolite in model.metabolites:
    globals()[metabolite.id] = metabolite

"""# Intracellular Core Reactions

TCA
"""

reaction = Reaction('Citrate_synthase')
reaction.name = 'Citrate synthase' #EC 2.3.3.1/2.3.3.3/2.3.3.16
reaction.subsystem = 'TCA'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({Acetyl_CoA: -1.0, Oxaloacetate: -1.0, Water_c: -1.0, Citrate: 1.0, CoA: 1.0, Proton_c: 1.0})
reaction.gene_reaction_rule = '(A8G16_RS23015 or A8G16_RS23035)'
model.add_reactions([reaction])

reaction = Reaction('Aconitase_A')
reaction.name = 'Aconitate hydratase A' #EC 4.2.1.3
reaction.subsystem = 'TCA'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({Citrate: -1.0, Cis_Aconitate: 1.0, Water_c: 1.0})
reaction.gene_reaction_rule = '(A8G16_RS08250 or A8G16_RS21735 or A8G16_RS39530)'
model.add_reactions([reaction])

reaction = Reaction('Aconitase_B') #EC 4.2.1.3
reaction.name = 'Aconitate hydratase B '
reaction.subsystem = 'TCA'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({Cis_Aconitate: -1.0, Water_c: -1.0, Isocitrate: 1.0})
reaction.gene_reaction_rule = '(A8G16_RS08250 or A8G16_RS21735 or A8G16_RS39530)'
model.add_reactions([reaction])

reaction = Reaction('Isocitrate_dehydrogenase_NADP')
reaction.name = 'Isocitrate dehydrogenase (NADP-dependent)' #EC 1.1.1.42
reaction.subsystem = 'TCA'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({Isocitrate: -1.0, NADP: -1.0, Two_Ketoglutarate: 1.0, CO2_c: 1.0, NADPH: 1.0})
reaction.gene_reaction_rule = '(A8G16_RS25530 or A8G16_RS28305)'
model.add_reactions([reaction])

reaction = Reaction('Two_ketoglutarate_dehydrogenase')
reaction.name = '2-ketoglutarate dehydrogenase' 
reaction.subsystem = 'TCA'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({Two_Ketoglutarate: -1.0, CoA: -1.0, NAD: -1.0, CO2_c: 1.0, Succinyl_CoA: 1.0, NADH: 1.0})
reaction.gene_reaction_rule = '((A8G16_RS01180 or A8G16_RS14355	or A8G16_RS21205	or A8G16_RS23340	or A8G16_RS28185) and A8G16_RS12715 and A8G16_RS24815)'
model.add_reactions([reaction])

reaction = Reaction('Succinyl_CoA_synthase') #EC 6.2.1.5
reaction.name = 'Succinyl-CoA synthetase'
reaction.subsystem = 'TCA'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({ADP: -1.0, Phosphate: -1.0, Succinyl_CoA: -1.0, ATP: 1.0, CoA: 1.0, Succinate: 1.0})
reaction.gene_reaction_rule = '(A8G16_RS32760 and A8G16_RS32765)'
model.add_reactions([reaction])

reaction = Reaction('Succinate_dehydrogenase')
reaction.name = 'Succinate dehydrogenase (ubiquinone)'
reaction.subsystem = 'TCA'
reaction.lower_bound = -1000. #EC 1.3.5.1
reaction.upper_bound = 1000.
reaction.add_metabolites({Succinate: -1.0, Ubiquinone8: -1.0, Ubiquinol8: 1.0, Fumarate: 1.0})
reaction.gene_reaction_rule = '(A8G16_RS11055 and A8G16_RS11060 and A8G16_RS11065 and A8G16_RS11070 )'
model.add_reactions([reaction])

reaction = Reaction('Fumarase')
reaction.name = 'Fumarase' #EC 4.2.1.2
reaction.subsystem = 'TCA'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({Fumarate: -1.0, Water_c: -1.0, Malate: 1.0})
reaction.gene_reaction_rule = '(A8G16_RS03235 or A8G16_RS03410)'
model.add_reactions([reaction])

reaction = Reaction('Malate_dehydrogenase_NAD')
reaction.name = 'Malate dehydrogenase NAD' #EC 1.1.1.37
reaction.subsystem = 'TCA'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({Malate: -1.0, NAD: -1.0, Oxaloacetate: 1.0, NADH: 1.0, Proton_c: 1.0})
reaction.gene_reaction_rule = '(A8G16_RS11075)'
model.add_reactions([reaction])


"""GLYOXYLATE SHUNT"""

reaction = Reaction('Isocitrate_lyase')
reaction.name = 'Isocitrate lyase' #EC 4.1.3.1
reaction.subsystem = 'GLYX Shunt'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({Isocitrate: -1.0, Succinate: 1.0, Glyoxylate: 1.0})
reaction.gene_reaction_rule = '(A8G16_RS01270)'
model.add_reactions([reaction])

reaction = Reaction('Malate_synthase')
reaction.name = 'Malate synthase' #EC 2.3.3.9
reaction.subsystem = 'GLYX Shunt'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({Acetyl_CoA: -1.0, Glyoxylate: -1.0, Water_c: -1.0, CoA: 1.0, Malate: 1.0, Proton_c: 1.0})
reaction.gene_reaction_rule = '(A8G16_RS38450 or A8G16_RS28620)'
model.add_reactions([reaction])

"""Anaplerotic reactions"""

reaction = Reaction('Phosphoenolpyruvate_carboxykinase_GTP') #EC  4.1.1.32 GTP
reaction.name = 'Phosphoenolpyruvate carboxykinase (GTP)'
reaction.subsystem = 'Anaplerotic Pathway'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({GTP: -1.0, Oxaloacetate: -1.0, GDP:1.0, CO2_c:1.0, Phosphoenolpyruvate:1.0})
reaction.gene_reaction_rule = '(A8G16_RS15300 or A8G16_RS29775 or A8G16_RS30415 or A8G16_RS31685)'
model.add_reactions([reaction])


reaction = Reaction('Pyruvate_carboxylase') #EC 6.4.1.1 
reaction.name = 'Pyruvate carboxylase	'
reaction.subsystem = 'Anaplerotic Pathway'
reaction.lower_bound = 0
reaction.upper_bound = 1000.
reaction.add_metabolites({ATP: -1.0, Pyruvate:-1.0, Hydrogen_carbonate: -1.0, Oxaloacetate: 1.0, ADP:1.0, Phosphate: 1.0, Proton_c:1.0})
reaction.gene_reaction_rule = '(A8G16_RS22465)'
model.add_reactions([reaction])

reaction = Reaction('Malic_Enzyme_NAD')
reaction.name = 'Malic enzyme NAD' #EC 1.1.1.38/1.1.1.39
reaction.subsystem = 'Anaplerotic Pathway'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({Malate: -1.0, NAD: -1.0, NADH:1.0, CO2_c:1.0, Pyruvate:1.0})
reaction.gene_reaction_rule = '(A8G16_RS16130 or A8G16_RS29645 or A8G16_RS32440)'
model.add_reactions([reaction])

reaction = Reaction('Malic_Enzyme_NADP')
reaction.name = 'Malic enzyme NADP' #EC 1.1.1.40
reaction.subsystem = 'Anaplerotic Pathway'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({Malate: -1.0, NADP: -1.0, NADPH:1.0, CO2_c:1.0, Pyruvate:1.0})
reaction.gene_reaction_rule = '(A8G16_RS14495 or A8G16_RS24755 or A8G16_RS28700 or A8G16_RS38355)'
model.add_reactions([reaction])

reaction = Reaction('Phosphoenolpyruvate_carboxylase')
reaction.name = 'Phosphoenolpyruvate carboxylase' #EC 4.1.1.31 
reaction.subsystem = 'Anaplerotic Pathway '
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({Phosphoenolpyruvate: -1.0, Hydrogen_carbonate: -1.0, Oxaloacetate:1.0, Phosphate:1.0 })
reaction.gene_reaction_rule = '(A8G16_RS33825)'
model.add_reactions([reaction])

"""Electron transport Chain, ATP synthesis and maintenance. And anothers rxns related to energy metabolism"""

reaction = Reaction('NADH_dehydrogenase')
reaction.name = 'NADH dehydrogenase'
reaction.subsystem = 'ETC' #EC 7.1.1.2
reaction.lower_bound = 0.0
reaction.upper_bound = 1000
reaction.add_metabolites({NADH: -1.0, Proton_c: -5.0, Ubiquinone8: -1.0, Ubiquinol8: 1.0, NAD: 1.0, Proton_e: 4.0})
reaction.gene_reaction_rule = '(A8G16_RS03185 and A8G16_RS03180 and A8G16_RS03175 and A8G16_RS03170 and A8G16_RS03165 and A8G16_RS03160 and A8G16_RS03155 and A8G16_RS03150 and A8G16_RS03145 and A8G16_RS03140 and A8G16_RS03135 and A8G16_RS03130 and A8G16_RS03125 and A8G16_RS03120)'
model.add_reactions([reaction])

reaction = Reaction('Oxygen_oxidoreductase')
reaction.name = 'Ubiquinol: oxygen oxidoreductase'
reaction.subsystem = 'ETC' #EC 7.1.1.3
reaction.lower_bound = 0.0
reaction.upper_bound = 1000.0
reaction.add_metabolites({Ubiquinol8: -2.0, Oxygen_c: -1.0, Proton_c: -5.0, Proton_e: 5.0, Ubiquinone8: 2.0, Water_c: 2.0})
reaction.gene_reaction_rule = '(A8G16_RS09830)'
model.add_reactions([reaction])

reaction = Reaction('ATP_synthase')
reaction.name = 'ATP synthase'
reaction.subsystem = 'ATP Biosynthesis'
reaction.lower_bound = 0
reaction.upper_bound = 1000.0
reaction.add_metabolites({ADP: -1.0, Proton_e: -5.0, Phosphate: -1.0, ATP: 1.0, Proton_c: 4.0, Water_c: 1.0})
reaction.gene_reaction_rule = '(A8G16_RS32060 and A8G16_RS32065 and A8G16_RS32070 and A8G16_RS32075 and A8G16_RS32080 and A8G16_RS32085 and A8G16_RS32090 and A8G16_RS32095)'
model.add_reactions([reaction])

reaction = Reaction('ATPase')
reaction.name = 'ATP for maintenance'
reaction.subsystem = 'ATPm'
reaction.lower_bound = 6
reaction.upper_bound = 1000
reaction.add_metabolites({ATP: -1.0, Water_c: -1.0, ADP: 1.0, Proton_c: 1.0, Phosphate: 1.0})
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])


reaction = Reaction('Transhydrogenase')
reaction.name = 'Proton-translocating NAD(P)+ transhydrogenase'
reaction.subsystem = 'Transhydrogenase'
reaction.lower_bound = 0 
reaction.upper_bound = 0
reaction.add_metabolites({NADPH: -1.0, NAD: -1.0, NADH: 1.0, NADP: 1.0})
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])


reaction = Reaction('Inorganic_pyrophosphatase')
reaction.name = 'Inorganic pyrophosphatase'
reaction.subsystem = 'Inorganic phosphate synthesis'
reaction.lower_bound = 0.0
reaction.upper_bound = 1000.0
reaction.add_metabolites({Diphosphate: -1.0, Water_c: -1.0, Phosphate: 2.0, Proton_c: 1.0})
reaction.gene_reaction_rule = '(A8G16_RS04705 or A8G16_RS12275)'
model.add_reactions([reaction])

reaction = Reaction('Adenylate_kinase')
reaction.name = 'Adenylate kinase'
reaction.subsystem = 'Adenosine ribonucleotides de novo biosynthesis'
reaction.lower_bound = -1000.0
reaction.upper_bound = 1000.0
reaction.add_metabolites({AMP: -1.0, ATP: -1.0, ADP: 2.0})
reaction.gene_reaction_rule = '(A8G16_RS00860 or A8G16_RS36010)'
model.add_reactions([reaction])

reaction = Reaction('Nucleoside_diphosphate_kinase') #EC  2.7.4.6
reaction.name = 'Nucleoside-diphosphate kinase'
reaction.subsystem = 'Guanosine ribonucleotides de novo biosynthesis'
reaction.lower_bound = 0
reaction.upper_bound = 1000.0
reaction.add_metabolites({GDP: -1.0, ATP: -1.0, ADP: 1.0, GTP:1.0})
reaction.gene_reaction_rule = '(A8G16_RS15380)'
model.add_reactions([reaction])

"""# Synthesis of precusors of TAGs"""

reaction = Reaction('Carbonic_anhydrase')
reaction.name = 'Carbonic anhydrase' #EC 4.2.1.1
reaction.subsystem = 'Synthesis of precursors of TAGs'
reaction.lower_bound = -1000.0
reaction.upper_bound = 1000.0
reaction.add_metabolites({Hydrogen_carbonate: -1.0, Proton_c: -1.0, CO2_c: 1.0, Water_c:1.0})
reaction.gene_reaction_rule = '(A8G16_RS04445 or A8G16_RS21825 or A8G16_RS22500 or A8G16_RS25040 or A8G16_RS29755 or A8G16_RS33255 or A8G16_RS35570)'
model.add_reactions([reaction])

"""Propionyl-CoA Synthesis"""

reaction = Reaction('Propionyl_CoA_Synthesis')
reaction.name = 'Propionyl-CoA Synthesis Lumped reaction'
reaction.subsystem = 'Synthesis of precursors of TAGs'
reaction.lower_bound = -1000. 
reaction.upper_bound = 1000.
reaction.add_metabolites({Succinyl_CoA: -1.0, ADP: -1.0, Proton_c: -1.0, Phosphate:-1.0, Propionyl_CoA: 1.0, ATP: 1.0, Hydrogen_carbonate: 1.0})
reaction.gene_reaction_rule = '((A8G16_RS27670 or A8G16_RS27675) and A8G16_RS12140 and (A8G16_RS00430 or A8G16_RS35710	or A8G16_RS37355	or A8G16_RS37580))'
model.add_reactions([reaction])

reaction = Reaction('Acetyl_CoA_carboxylase')
reaction.name = 'Acetyl-CoA carboxylase' #EC 6.3.4.14, 2.1.3.15
reaction.subsystem = 'Synthesis of precursors of TAGs'
reaction.subsystem = 'Fatt y acid biosynthesis initation'
reaction.lower_bound = 0.0
reaction.upper_bound = 1000.0
reaction.add_metabolites({ATP: -1.0, Acetyl_CoA: -1.0, Hydrogen_carbonate: -1.0, ADP: 1.0, Malonyl_CoA: 1.0, Phosphate: 1.0, Proton_c: 1.0})
reaction.gene_reaction_rule = '(A8G16_RS18105 and A8G16_RS00425 and (A8G16_RS21065 or A8G16_RS25685))'
model.add_reactions([reaction])

reaction = Reaction('Two_oxoglutarate_synthase')
reaction.name = '	2-oxoglutarate synthase'
reaction.subsystem = 'Synthesis of precursors of TAGs' #EC 1.2.7.3 , ferredoxins synthesis
reaction.lower_bound = -1000.0
reaction.upper_bound = 1000.0
reaction.add_metabolites({CoA: -1.0, Oxidized_ferredoxins: -2.0, Two_Ketoglutarate: -1.0, Succinyl_CoA: 1.0, CO2_c: 1.0, Proton_c: 1.0, Reduced_ferredoxins: 2.0})
reaction.gene_reaction_rule = '(A8G16_RS08295 and A8G16_RS08300)'
model.add_reactions([reaction])

reaction = Reaction('C15_1_synthesis')
reaction.name = 'Pentadecenoic acid synthesis'
reaction.subsystem = 'Synthesis of precursors of TAGs'
reaction.lower_bound = 0.0
reaction.upper_bound = 1000.0
reaction.add_metabolites({Pentadecanoate: -1.0, Oxygen_c: -1.0, Proton_c: -2.0, Reduced_ferredoxins: -2.0, Pentadecenoic_acid: 1.0, Oxidized_ferredoxins: 2.0, Water_c: 2.0})
reaction.gene_reaction_rule = '(A8G16_RS30250 or A8G16_RS20145 or A8G16_RS20140 or A8G16_RS14585 or A8G16_RS12395 or A8G16_RS09935 or A8G16_RS09280 or A8G16_RS09095 or A8G16_RS08690 or A8G16_RS07190 or A8G16_RS04365 or A8G16_RS03415)'
model.add_reactions([reaction])

reaction = Reaction('C14_synthesis')
reaction.name = 'Tetradecanoic acid synthesis'
reaction.subsystem = 'Synthesis of precursors of TAGs'
reaction.lower_bound = 0.0
reaction.upper_bound = 1000.0
reaction.add_metabolites({Acetyl_CoA: -1.0, Malonyl_CoA: -6.0, NADPH: -12.0, Proton_c: -17.0, Tetradecanoate: 1.0, CoA: 7.0, CO2_c: 6.0, Water_c: 5.0, NADP: 12.0})
reaction.gene_reaction_rule = '(A8G16_RS12305)'
model.add_reactions([reaction])

reaction = Reaction('C16_synthesis')
reaction.name = 'Palmitic acid synthesis'
reaction.subsystem = 'Synthesis of precursors of TAGs'
reaction.lower_bound = 0.0
reaction.upper_bound = 1000.0
reaction.add_metabolites({Acetyl_CoA: -1.0, Malonyl_CoA: -7.0, NADPH: -14.0, Proton_c: -20.0, Palmitate: 1.0, CoA: 8.0, CO2_c: 7.0, Water_c: 6.0, NADP: 14.0})
reaction.gene_reaction_rule = '(A8G16_RS12305)'
model.add_reactions([reaction])

reaction = Reaction('C18_synthesis')
reaction.name = 'Stearic Acid synthesis'
reaction.subsystem = 'Synthesis of precursors of TAGs'
reaction.lower_bound = 0.0
reaction.upper_bound = 1000.0
reaction.add_metabolites({Acetyl_CoA: -1.0, Malonyl_CoA: -8.0, NADPH: -16.0, Proton_c: -23.0, Stereate: 1.0, CoA: 9.0, CO2_c: 8.0, Water_c: 7.0, NADP: 16.0})
reaction.gene_reaction_rule = '(A8G16_RS12305)'
model.add_reactions([reaction])

reaction = Reaction('C15_sythesis')
reaction.name = 'Pentadecanoic acid synthesis'
reaction.subsystem = 'Synthesis of precursors of TAGs'
reaction.lower_bound = 0.0
reaction.upper_bound = 1000.0
reaction.add_metabolites({Propionyl_CoA: -1.0, Malonyl_CoA: -6.0, NADPH: -12.0, Proton_c: -17.0, Pentadecanoate: 1.0, CoA: 7.0, CO2_c: 6.0, Water_c: 5.0, NADP: 12.0})
reaction.gene_reaction_rule = '(A8G16_RS12305)'
model.add_reactions([reaction])

reaction = Reaction('C17_synthesis')
reaction.name = 'Heptadecanoic acid/Margaric acid synthesis'
reaction.subsystem = 'Synthesis of precursors of TAGs'
reaction.lower_bound = 0.0
reaction.upper_bound = 1000.0
reaction.add_metabolites({Pentadecanoate: -1.0, Malonyl_CoA: -1.0, NADPH: -2.0, Proton_c: -3.0, Heptadecanoate: 1.0, CoA: 1.0, CO2_c: 1.0, Water_c: 1.0, NADP: 2.0})
reaction.gene_reaction_rule = '(A8G16_RS12305)'
model.add_reactions([reaction])

reaction = Reaction('C17_1_synthesis')
reaction.name = 'Heptadecenoic acid syhtnesis'
reaction.subsystem = 'Synthesis of precursors of TAGs'
reaction.lower_bound = 0.0
reaction.upper_bound = 1000.0
reaction.add_metabolites({Heptadecanoate: -1.0, Oxygen_c: -1.0, Proton_c: -2.0, Reduced_ferredoxins: -2.0, Margaroleate: 1.0, Oxidized_ferredoxins: 2.0, Water_c: 2.0})
reaction.gene_reaction_rule = '(A8G16_RS30250 or A8G16_RS20145 or A8G16_RS20140 or A8G16_RS14585 or A8G16_RS12395 or A8G16_RS09935 or A8G16_RS09280 or A8G16_RS09095 or A8G16_RS08690 or A8G16_RS07190 or A8G16_RS04365 or A8G16_RS03415)'
model.add_reactions([reaction])

reaction = Reaction('C16_1_synthesis')
reaction.name = 'Palmitoleic acid synthesis'
reaction.subsystem = 'Synthesis of precursors of TAGs'
reaction.lower_bound = 0.0
reaction.upper_bound = 1000.0
reaction.add_metabolites({Palmitate: -1.0, Oxygen_c: -1.0, Proton_c: -2.0, Reduced_ferredoxins: -2.0, Palmitoleate: 1.0, Oxidized_ferredoxins: 2.0, Water_c: 2.0})
reaction.gene_reaction_rule = '(A8G16_RS30250 or A8G16_RS20145 or A8G16_RS20140 or A8G16_RS14585 or A8G16_RS12395 or A8G16_RS09935 or A8G16_RS09280 or A8G16_RS09095 or A8G16_RS08690 or A8G16_RS07190 or A8G16_RS04365 or A8G16_RS03415)'
model.add_reactions([reaction])

reaction = Reaction('C18_1_synthesis')
reaction.name = 'Oleic acid synthesis'
reaction.subsystem = 'Synthesis of precursors of TAGs'
reaction.lower_bound = 0.0
reaction.upper_bound = 1000.0
reaction.add_metabolites({Stereate: -1.0, Oxygen_c: -1.0, Proton_c: -2.0, Reduced_ferredoxins: -2.0, Cis_Trans_Vaccenate: 1.0, Oxidized_ferredoxins: 2.0, Water_c: 2.0})
reaction.gene_reaction_rule = '(A8G16_RS30250 or A8G16_RS20145 or A8G16_RS20140 or A8G16_RS14585 or A8G16_RS12395 or A8G16_RS09935 or A8G16_RS09280 or A8G16_RS09095 or A8G16_RS08690 or A8G16_RS07190 or A8G16_RS04365 or A8G16_RS03415)'
model.add_reactions([reaction])

"""# Assimilation of substrates

Glucose
"""

#IMPORT OF GLUCOSE WITH PTS
reaction = Reaction('PTS_glucose_transport')
reaction.name = 'PTS glucose/sucrose transporter'
reaction.subsystem = 'EMP/Glucose transport'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({Glucose_e: -1.0, Phosphoenolpyruvate: -1.0, Pyruvate: 1.0, D_Glucose_6_phosphate:1.0})
reaction.gene_reaction_rule = '(A8G16_RS29820 and A8G16_RS29825 and A8G16_RS29830 and A8G16_RS29835 and A8G16_RS29840)'
model.add_reactions([reaction])

reaction = Reaction('Glucose_6_phosphate_isomerase')
reaction.name = 'Glucose-6-phosphate isomerase'
reaction.subsystem = 'EMP'
reaction.lower_bound = -1000 # EC 5.3.1.9
reaction.upper_bound = 1000
reaction.add_metabolites({D_Glucose_6_phosphate: -1.0, D_Fructose_6_phosphate: 1.0})
reaction.gene_reaction_rule = '(A8G16_RS00040 or A8G16_RS39595)'
model.add_reactions([reaction])

reaction = Reaction('six_phosphofructokinase')
reaction.name = '6-phosphofructokinase'
reaction.subsystem = 'EMP'
reaction.lower_bound = 0. # EC 2.7.1.11 
reaction.upper_bound = 1000
reaction.add_metabolites({D_Fructose_6_phosphate: -1.0, ATP: -1.0, ADP:1.0, Proton_c:1.0, D_Fructose_1_6_bisphosphate:1.0})
reaction.gene_reaction_rule = '(A8G16_RS11775)'
model.add_reactions([reaction])

reaction = Reaction('Fructose_bisphosphate_aldolase')
reaction.name = 'Fructose bisphosphate aldolase'
reaction.subsystem = 'EMP/CBB' #EC 4.1.2.13
reaction.lower_bound = 0.#-1000
reaction.upper_bound = 1000.
reaction.add_metabolites({D_Fructose_1_6_bisphosphate: -1.0, Glyceraldehyde_3_phosphate: 1.0, Dihydroxyacetone_phosphate:1.0})
reaction.gene_reaction_rule = '(A8G16_RS25905 or A8G16_RS30370 or A8G16_RS41170)'
model.add_reactions([reaction])

reaction = Reaction('Triose_phosphate_isomerase')
reaction.name = 'Triose-phosphate isomerase'
reaction.subsystem = 'EMP/CBB'
reaction.lower_bound = -1000. # EC 5.3.1.1
reaction.upper_bound = 1000.
reaction.add_metabolites({Glyceraldehyde_3_phosphate: 1.0, Dihydroxyacetone_phosphate:-1.0})
reaction.gene_reaction_rule = '(A8G16_RS33835)'
model.add_reactions([reaction])


reaction = Reaction('glyceraldehyde_3_phosphate_dehydrogenase_NADP')
reaction.name = 'Glyceraldehyde-3-phosphate dehydrogenase (NADP+) (phosphorylating)'
reaction.subsystem = 'EMP/CBB' #GAPDHN 1.2.1.13 (NADP)
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({Glyceraldehyde_3_phosphate: -1.0, NADP:-1.0, Phosphate:-1.0, Proton_c:1.0, NADPH:1.0, Phospho3_D_glyceroyl_phosphate:1.0})
reaction.gene_reaction_rule = '(A8G16_RS02100)'
model.add_reactions([reaction])

reaction = Reaction('Phosphoglycerate_kinase')
reaction.name = 'phosphoglycerate kinase'
reaction.subsystem = 'EMP/CBB' # EC 2.7.2.3
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({Phospho3_D_glycerate:1.0, ATP:1.0, ADP:-1.0, Phospho3_D_glyceroyl_phosphate:-1.0})
reaction.gene_reaction_rule = '(A8G16_RS25850 or A8G16_RS33840)'
model.add_reactions([reaction])

reaction = Reaction('Phosphoglycerate_mutase')
reaction.name = 'phosphoglycerate mutase'
reaction.subsystem = 'EMP' #EC 5.4.2.11
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({Phospho3_D_glycerate:-1.0, D_Glycerate_2_phosphate:1.0})
reaction.gene_reaction_rule = '(A8G16_RS01505)'
model.add_reactions([reaction])

reaction = Reaction('Enolase')
reaction.name = 'Enolase'
reaction.subsystem = 'EMP' #EC 4.2.1.11
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({Phosphoenolpyruvate:1.0, Water_c: 1.0,D_Glycerate_2_phosphate:-1.0})
reaction.gene_reaction_rule = '(A8G16_RS14105 or A8G16_RS24170)'
model.add_reactions([reaction])

reaction = Reaction('Pyruvate_kinase')
reaction.name = 'Pyruvate kinase'
reaction.subsystem = 'EMP' #EC 2.7.1.40
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({Phosphoenolpyruvate:-1.0, ADP: -1.0,Proton_c:-1.0, ATP:1.0, Pyruvate:1.0})
reaction.gene_reaction_rule = '(A8G16_RS09620 or A8G16_RS14250 or A8G16_RS27205 or A8G16_RS39215 or A8G16_RS42505)'
model.add_reactions([reaction])

reaction = Reaction('Glycerol_3_phosphate_dehydrogenase_NADP')
reaction.name = 'Glycerol-3-phosphate dehydrogenase NADP'
reaction.subsystem = 'Synthesis of precursors of TAGs' #EC 1.1.1.94
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({Dihydroxyacetone_phosphate:-1.0, NADPH: -1.0,Proton_c:-1.0, NADP:1.0, Glycerol_3P:1.0})
reaction.gene_reaction_rule = '(A8G16_RS08265 or A8G16_RS33225)'
model.add_reactions([reaction])

reaction = Reaction('Pyruvate_dehydrogenase_NAD')
reaction.name = 'pyruvate dehydrogenase NAD'
reaction.subsystem = 'Pyruvate decarboxylation to acetyl CoA' #EC 1.2.1.104
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({Pyruvate:-1.0, NAD: -1.0,CoA:-1.0, NADH:1.0, Acetyl_CoA:1.0, CO2_c:1.0})
reaction.gene_reaction_rule = '((A8G16_RS01180 or A8G16_RS14355	or A8G16_RS21205	or A8G16_RS23340	or A8G16_RS28185) and (A8G16_RS04155 or A8G16_RS36855) and (A8G16_RS07575 or A8G16_RS10630 or A8G16_RS14370 or A8G16_RS28170 or A8G16_RS04155))'
model.add_reactions([reaction])

reaction = Reaction('Phosphoenolpyruvate_synthase')
reaction.name = 'Phosphoenolpyruvate_synthase'
reaction.subsystem = 'Gluconeogenesis' # EC 2.7.9.2
reaction.lower_bound = 0
reaction.upper_bound = 1000.0
reaction.add_metabolites({ATP: -1.0, Pyruvate: -1.0, Water_c: -1.0, AMP: 1.0, Phosphate: 1.0, Phosphoenolpyruvate: 1.0, Proton_c: 2.0})
reaction.gene_reaction_rule = '(A8G16_RS10130 or A8G16_RS07615)'
model.add_reactions([reaction])

"""EPS synthesis"""

reaction = Reaction('Glucose_phosphomutase')
reaction.name = 'Glucose phosphomutase'
reaction.subsystem = 'Biosynthesis EPS' #EC 5.4.2.2
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({Glucose_1_phosphate: 1.0, D_Glucose_6_phosphate:-1.0 })
reaction.gene_reaction_rule = '(A8G16_RS21880)'
model.add_reactions([reaction])

reaction = Reaction('Glucose_1_phosphate_adenylyltransferase')
reaction.name = 'Glucose 1-phosphate adenylyltransferase'
reaction.subsystem = 'Biosynthesis EPS' #EC 2.7.7.27
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({Glucose_1_phosphate: -1.0, ATP:-1.0, Proton_c:-1.0, ADPglucose:1.0, Diphosphate:1.0})
reaction.gene_reaction_rule = '(A8G16_RS24285)'
model.add_reactions([reaction])

reaction = Reaction('Glycogen_synthase')
reaction.name = 'Glycogen synthase'
reaction.subsystem = 'Biosynthesis EPS'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({ADPglucose: -1.0, ADP:1.0, Proton_c:1.0, Glycogen_c:1.0})
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

"""ACETATE CONVERSION TO ACETYL-COA"""

reaction = Reaction('Acetate_kinase')
reaction.name = 'Acetate kinase'
reaction.subsystem = 'Acetate assimilation' #EC 2.7.2.1/2.7.2.15
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({Acetate_c: -1.0, ATP: -1.0, ADP: 1.0, Acetyl_phosphate: 1.0})
reaction.gene_reaction_rule = '(A8G16_RS34015)'
model.add_reactions([reaction])

reaction = Reaction('Phosphate_acetyltransferase')
reaction.name = 'Phosphate acetyltransferase'
reaction.subsystem = 'Acetate assimilation' #EC 2.3.1.8
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({Acetyl_phosphate: -1.0, CoA: -1.0, Phosphate:1.0, Acetyl_CoA: 1.0})
reaction.gene_reaction_rule = '(A8G16_RS34010)'
model.add_reactions([reaction])

"""CALVIN-BENSON-BASSHAM"""

reaction = Reaction('RuBisCO')
reaction.name = 'Ribulose bisphosphate carboxylase' #EC 4.1.1.39
reaction.subsystem = 'CBB'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({D_Ribulose_1_5_bisphosphate: -1.0, CO2_c: -1.0, Water_c:-1.0, Phospho3_D_glycerate:2.0, Proton_c:2.0})
reaction.gene_reaction_rule = '(A8G16_RS25910 and A8G16_RS25915)'
model.add_reactions([reaction])

reaction = Reaction('Fructose_1_6_diphosphatase')
reaction.name = 'fructose 1,6-diphosphatase'
reaction.subsystem = 'CBB/Gluconeogenesis' # EC 3.1.3.11
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({D_Fructose_1_6_bisphosphate: -1.0, Water_c: -1.0, Phosphate:1.0, D_Fructose_6_phosphate:1.0})
reaction.gene_reaction_rule = '(A8G16_RS00015 or A8G16_RS03405 or A8G16_RS25835)'
model.add_reactions([reaction])

reaction = Reaction('TKT2')
reaction.name = 'transketolase 2'
reaction.subsystem = 'CBB/PPP'
reaction.lower_bound = -1000. # EC 2.2.1.1
reaction.upper_bound = 1000.
reaction.add_metabolites({D_Fructose_6_phosphate:-1.0, Glyceraldehyde_3_phosphate:-1.0, D_Erythrose_4_phosphate: 1.0, D_Xylulose_5_phosphate:1.0 })
reaction.gene_reaction_rule = '(A8G16_RS25895)'
model.add_reactions([reaction])

reaction = Reaction('Sedoheptulose_1_7_bisphosphate_aldolase')
reaction.name = 'Fructose-bisphosphate aldolase/sedoheptulose-1,7-bisphosphate aldolase'
reaction.subsystem = 'CBB/PPP' # EC 4.1.2.13
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.
reaction.add_metabolites({D_Erythrose_4_phosphate: -1.0, Dihydroxyacetone_phosphate: -1.0, Sedoheptulose_1_7_bisphosphate:1.0})
reaction.gene_reaction_rule = '(A8G16_RS25905 or A8G16_RS30370 or A8G16_RS41170)'
model.add_reactions([reaction])

reaction = Reaction('Sedoheptulose_bisphosphatase')
reaction.name = 'Sedoheptulose-bisphosphatase'
reaction.subsystem = 'CBB'
reaction.lower_bound = 0. # EC 3.1.3.11
reaction.upper_bound = 1000.
reaction.add_metabolites({Sedoheptulose_1_7_bisphosphate:-1.0, Water_c:-1.0, Phosphate:1.0, Sedoheptulose_7_phosphate:1.0 })
reaction.gene_reaction_rule = '(A8G16_RS00015 or A8G16_RS03405 or A8G16_RS25835)'
model.add_reactions([reaction])

reaction = Reaction('TKT1')
reaction.name = 'transketolase 1'
reaction.subsystem = 'CBB/PPP'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000. # EC 2.2.1.1
reaction.add_metabolites({Sedoheptulose_7_phosphate:-1.0, Glyceraldehyde_3_phosphate:-1.0, D_Xylulose_5_phosphate:1.0, Alpha_D_Ribose_5_phosphate:1.0 })
reaction.gene_reaction_rule = '(A8G16_RS25895)'
model.add_reactions([reaction])

reaction = Reaction('Ribose_5_phosphate_isomerase')
reaction.name = 'Ribose-5-phosphate isomerase'
reaction.subsystem = 'CBB/PPP'
reaction.lower_bound = -1000. #EC 5.3.1.6
reaction.upper_bound = 1000.
reaction.add_metabolites({Alpha_D_Ribose_5_phosphate:-1.0, D_Ribulose_5_phosphate:1.0 })
reaction.gene_reaction_rule = '(A8G16_RS23885)'
model.add_reactions([reaction])

reaction = Reaction('Phosphoribulokinase')
reaction.name = 'Phosphoribulokinase'
reaction.subsystem = 'CBB'
reaction.lower_bound = 0.
reaction.upper_bound = 1000. # EC 2.7.1.19
reaction.add_metabolites({D_Ribulose_5_phosphate:-1.0, ATP:-1.0, ADP:1.0, Proton_c:1.0, D_Ribulose_1_5_bisphosphate:1.0 })
reaction.gene_reaction_rule = '(A8G16_RS25840)'
model.add_reactions([reaction])

reaction = Reaction('Ribulose_phosphate_3_epimerase')
reaction.name = 'Ribulose phosphate 3-epimerase'
reaction.subsystem = 'CBB/PPP'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000. # EC 5.1.3.1
reaction.add_metabolites({D_Ribulose_5_phosphate:-1.0, D_Xylulose_5_phosphate:1.0 })
reaction.gene_reaction_rule = '(A8G16_RS15870 or A8G16_RS25890)'
model.add_reactions([reaction])

"""Hydrogenase

"""

reaction = Reaction('Hydrogenase_NADH_dependent')
reaction.name = 'NADH-dependent hydrogenase'
reaction.subsystem = 'Hydrogen assimilation' #EC 1.12.1.2
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({NAD:-1.0, H2_c:-1.0, Proton_c:1.0, NADH:1.0})
reaction.gene_reaction_rule = '(A8G16_RS25965 and A8G16_RS25960 and A8G16_RS25955 and A8G16_RS25950 and A8G16_RS25945 and A8G16_RS25940)'
model.add_reactions([reaction])

"""Propanol assimilation"""

reaction = Reaction('Propanol_Assimilation')
reaction.name = 'Lumped reaction: 1-Propanol Assimilation'
reaction.subsystem = '1-Propanol Assimilation'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({NAD:-1.0, Propanol_c:-1.0, CoA:-1.0, NADP:-1.0, ATP:-1.0, Water_c:-1.0, Propionyl_CoA:1.0, NADH:1.0, NADPH:1.0, AMP:1.0, Diphosphate: 1.0, Proton_c: 3.0})
reaction.gene_reaction_rule = '((A8G16_RS09425 or A8G16_RS09950	 or A8G16_RS23305	 or A8G16_RS24205	 or A8G16_RS30960	 or A8G16_RS34995	or A8G16_RS36340 or A8G16_RS36405) and (A8G16_RS08715	or A8G16_RS14530	or A8G16_RS32285	or A8G16_RS34295) and (A8G16_RS02080	or A8G16_RS18080))'
model.add_reactions([reaction])

#Butanoate Assimilation
reaction = Reaction('Butanoyl_CoA_synthase')
reaction.name = 'Butanoyl_CoA_synthase' #EC 6.2.1.2
reaction.subsystem = 'Butanoate Assimilation'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({Butanoate_c:-1.0, CoA:-1.0, ATP:-1.0, Butanoyl_CoA:1.0, AMP:1.0, Diphosphate:1.0})
reaction.gene_reaction_rule = '(A8G16_RS09430 or A8G16_RS19465 or A8G16_RS29970)'
model.add_reactions([reaction])

reaction = Reaction('Butanoyl_CoA_assimilation')
reaction.name = 'Lumped reaction: Butanoyl-CoA to Acetyl-CoA'
reaction.subsystem = 'Butanoate Assimilation'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({Butanoyl_CoA:-1.0, NAD:-2.0, Water_c:-1.0, CoA:-1.0, Acetyl_CoA:2.0, NADH:2.0, Proton_c:2.0})
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

#Valerate

reaction = Reaction('Valerly_CoA_synthase')
reaction.name = 'Valerly CoA synthase'
reaction.subsystem = 'Valertate Assimilation'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({Valerate_c:-1.0, CoA:-1.0, ATP:-1.0, Valerly_CoA:1.0, AMP:1.0, Diphosphate:1.0})
reaction.gene_reaction_rule = '(A8G16_RS09430 or A8G16_RS19465 or A8G16_RS29970)'
model.add_reactions([reaction])


reaction = Reaction('Valeryl_CoA_assimilation')
reaction.name = 'Lumped reaction: Valeryl-CoA to Acetyl-CoA and Propionyl-CoA'
reaction.subsystem = 'Valertate Assimilation'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({Valerly_CoA:-1.0, NAD:-2.0, Water_c:-1.0, CoA:-1.0, Acetyl_CoA:1.0, NADH:2.0, Proton_c:2.0, Propionyl_CoA:1.0})
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

#Methanol

reaction = Reaction('Methanol_dehydrogenase_NAD')
reaction.name = 'Methanol dehydogenase NAD'
reaction.subsystem = 'Methanol assimilation to formaldehyde NAD+'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({Methanol_c:-1.0, NAD:-1.0, NADH:1.0, Formaldehyde:1.0, Proton_c:1.0})
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

reaction = Reaction('Formaldehyde_RuMP')
reaction.name = 'Lumped reaction: Formaldelhyde assimilation RuMP'
reaction.subsystem = 'Formaldehyde assimilation via RuMP' #Lumped
reaction.lower_bound = 0.
reaction.upper_bound = 0.
reaction.add_metabolites({Formaldehyde:-1.0, D_Ribulose_5_phosphate:-1.0, D_Fructose_6_phosphate:1.0 })
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

reaction = Reaction('Formaldehyde_Serine') #Lumped
reaction.name = 'Lumped reaction: Formaldelhyde assimilation Via Serine'
reaction.subsystem = 'Formaldehyde assimilation via Serine Cycle (THF)'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({Formaldehyde:-1.0, Glyoxylate:-1.0,NADPH:-1.0, ATP:-1.0, D_Glycerate_2_phosphate:1.0, ADP:1.0, NADP:1.0 })
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

reaction = Reaction('Formaldeyde_dissasimilation')
reaction.name = 'Lumped reaction: Formaldehyde, Formate dehydrogenase NAD' #Lumped
reaction.subsystem = 'Formaldehyde Dissimilation to CO2 '
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({Formaldehyde:-1.0, NAD:-2.0,Water_c:-1.0, NADH:2.0, CO2_c:1.0, Proton_c:2.0 })
reaction.gene_reaction_rule = '((A8G16_RS35155) and (A8G16_RS16200 and A8G16_RS16195))'
model.add_reactions([reaction])

reaction = Reaction('Glucose_6_phosphate_dehydrogenase_NADP')
reaction.name = 'NADP+-dependent glucose-6-phosphate dehydrogenase'
reaction.subsystem = 'ED' #EC 1.1.1.49 
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({D_Glucose_6_phosphate: -1.0, NADP: -1.0, Proton_c: 1.0, NADPH:1.0, D_glucono_1_5_lactone_6_phosphate:1.0})
reaction.gene_reaction_rule = '(A8G16_RS33805 or A8G16_RS13485 or A8G16_RS08400)'
model.add_reactions([reaction])

reaction = Reaction('Six-phosphogluconolactonase')
reaction.name = '6-phosphogluconolactonase' #EC 3.1.1.31 
reaction.subsystem = 'ED'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({D_glucono_1_5_lactone_6_phosphate: -1.0, Water_c: -1.0, Proton_c: 1.0, D_gluconate_6_phosphate:1.0})
reaction.gene_reaction_rule = '(A8G16_RS33815 or A8G16_RS00100)'
model.add_reactions([reaction])


reaction = Reaction('Phosphogluconate_dehydratase')
reaction.name = 'phosphogluconate dehydratase' # EC 4.2.1.12
reaction.subsystem = 'ED'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({D_gluconate_6_phosphate: -1.0, Water_c: 1.0, KDPG:1.0 })
reaction.gene_reaction_rule = '(A8G16_RS08405)'
model.add_reactions([reaction])

reaction = Reaction('KDPG_aldolase')
reaction.name = '2-dehydro-3-deoxy-phosphogluconate aldolase' # EC 4.1.2.14/4.1.2.55
reaction.subsystem = 'ED'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.
reaction.add_metabolites({KDPG: -1.0, Pyruvate: 1.0, Glyceraldehyde_3_phosphate:1.0 })
reaction.gene_reaction_rule = '(A8G16_RS08410)'
model.add_reactions([reaction])

reaction = Reaction('Transaldolase')
reaction.name = 'Transaldolase '
reaction.subsystem = 'PPP'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000. # EC 2.2.1.2
reaction.add_metabolites({D_Fructose_6_phosphate:-1.0, D_Erythrose_4_phosphate:-1.0, Sedoheptulose_7_phosphate:1.0, Glyceraldehyde_3_phosphate: 1.0})
reaction.gene_reaction_rule = '(A8G16_RS33800 or A8G16_RS21695 or A8G16_RS00035)'
model.add_reactions([reaction])

reaction = Reaction('Six_phosphogluconate_dehydrogenase')
reaction.name = '6-phosphogluconate dehydrogenase'
reaction.subsystem = 'PPP'
reaction.lower_bound = 0.
reaction.upper_bound = 1000. # 1.1.1.44
reaction.add_metabolites({D_gluconate_6_phosphate:-1.0, NADP:-1.0, D_Ribulose_5_phosphate:1.0, CO2_c: 1.0, NADPH:1.0})
reaction.gene_reaction_rule = '(A8G16_RS27605)'
model.add_reactions([reaction])

"""# TAGS reactions

Change the boundries depeding on the substrate
"""

reaction = Reaction('TAG-acetate')  
reaction.name = 'TAGs synthesis on acetate'
reaction.subsystem = 'Lipid Metabolism'
reaction.lower_bound = 0.0
reaction.upper_bound = 1000  # Irreversible

reaction.add_metabolites({
    ATP: -3.0,
    Margaroleate: -0.390,       # C17:1
    Heptadecanoate: -0.353,     # C17:0
    Tetradecanoate: -0.117,     # C14:0
    Pentadecanoate: -0.302,     # C15:0
    Pentadecenoic_acid: -0.024, # C15:1
    Palmitoleate: -0.301,       # C16:1
    Cis_Trans_Vaccenate: -0.411,# C18:1 (cis + trans)
    Glycerol_3P: -1.0,          # sn-Glycerol 3-phosphate
    Palmitate: -0.978,          # C16:0
    Stereate: -0.124,           # C18:0
    Water_c: -3.0,
    AMP: 3.0,
    Diphosphate: 3.0,
    Proton_c: 3.0,
    Phosphate: 1.0
})

reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

reaction = Reaction('TAG-1-Propanol')  
reaction.name = 'TAGs synthesis on 1-Propanol'
reaction.subsystem = 'Lipid Metabolism'
reaction.lower_bound = 0.0
reaction.upper_bound = 1000  # Irreversible

reaction.add_metabolites({
    ATP: -3.0,
    Margaroleate: -0.919,       # C17:1
    Heptadecanoate: -0.445,     # C17:0
    Tetradecanoate: -0.019,     # C14:0
    Pentadecanoate: -1.246,     # C15:0
    Pentadecenoic_acid: -0.161, # C15:1
    Palmitoleate: -0.055,       # C16:1
    Cis_Trans_Vaccenate: -0.054,# C18:1 (cis + trans)
    Glycerol_3P: -1.0,          # sn-Glycerol 3-phosphate
    Palmitate: -0.101,          # C16:0
    Stereate: -0.000,           # C18:0
    Water_c: -3.0,
    AMP: 3.0,
    Diphosphate: 3.0,
    Proton_c: 3.0,
    Phosphate: 1.0
})

reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

reaction = Reaction('TAG-ButiricAcid')
reaction.name = 'TAGs synthesis on Butiric Acid'
reaction.subsystem = 'Lipid Metabolism'
reaction.lower_bound = 0.0
reaction.upper_bound = 1000  # Irreversible

reaction.add_metabolites({
    ATP: -3.0,
    Margaroleate: -0.654,       # C17:1
    Heptadecanoate: -0.343,     # C17:0
    Tetradecanoate: -0.067,     # C14:0
    Pentadecanoate: -0.356,     # C15:0
    Pentadecenoic_acid: -0.028, # C15:1
    Palmitoleate: -0.334,       # C16:1
    Cis_Trans_Vaccenate: -0.458,# C18:1 (cis + trans)
    Glycerol_3P: -1.0,          # sn-Glycerol 3-phosphate
    Palmitate: -0.759,          # C16:0
    Stereate: -0.000,           # C18:0
    Water_c: -3.0,
    AMP: 3.0,
    Diphosphate: 3.0,
    Proton_c: 3.0,
    Phosphate: 1.0
})

reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

reaction = Reaction('TAG-ValericAcid')  
reaction.name = 'TAGs synthesis on Valeric Acid'
reaction.subsystem = 'Lipid Metabolism'
reaction.lower_bound = 0.0
reaction.upper_bound = 1000  # Irreversible

reaction.add_metabolites({
    ATP: -3.0,
    Margaroleate: -1.058,       # C17:1
    Heptadecanoate: -0.385,     # C17:0
    Tetradecanoate: -0.027,     # C14:0
    Pentadecanoate: -1.048,     # C15:0
    Pentadecenoic_acid: -0.152, # C15:1
    Palmitoleate: -0.095,       # C16:1
    Cis_Trans_Vaccenate: -0.084,# C18:1 (cis + trans)
    Glycerol_3P: -1.0,          # sn-Glycerol 3-phosphate
    Palmitate: -0.152,          # C16:0
    Stereate: -0.000,           # C18:0
    Water_c: -3.0,
    AMP: 3.0,
    Diphosphate: 3.0,
    Proton_c: 3.0,
    Phosphate: 1.0
})

reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

reaction = Reaction('TAG-Glucose')  
reaction.name = 'TAGs synthesis on Glucose'
reaction.subsystem = 'Lipid Metabolism'
reaction.lower_bound = 0
reaction.upper_bound = 1000  # Irreversible

reaction.add_metabolites({
    ATP: -3.0,
    Margaroleate: -0.482,       # C17:1
    Heptadecanoate: -0.339,     # C17:0
    Tetradecanoate: -0.116,     # C14:0
    Pentadecanoate: -0.235,     # C15:0
    Pentadecenoic_acid: -0.018, # C15:1
    Palmitoleate: -0.389,       # C16:1
    Cis_Trans_Vaccenate: -0.000,# C18:1 (cis + trans)
    Glycerol_3P: -1.0,          # sn-Glycerol 3-phosphate
    Palmitate: -1.240,          # C16:0
    Stereate: -0.181,           # C18:0
    Water_c: -3.0,
    AMP: 3.0,
    Diphosphate: 3.0,
    Proton_c: 3.0,
    Phosphate: 1.0
})

reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

reaction = Reaction('TAG-Knallgas')  
reaction.name = 'TAGs synthesis on CO2/H2/O2'
reaction.subsystem = 'Lipid Metabolism'
reaction.lower_bound = 0.0
reaction.upper_bound = 1000  # Irreversible

reaction.add_metabolites({
    ATP: -3.0,
    Margaroleate: -0.375,       # C17:1
    Heptadecanoate: -0.216,     # C17:0
    Tetradecanoate: -0.127,     # C14:0
    Pentadecanoate: -0.272,     # C15:0
    Pentadecenoic_acid: -0.000, # C15:1
    Palmitoleate: -0.484,       # C16:1
    Cis_Trans_Vaccenate: -0.435,# C18:1 (cis + trans)
    Glycerol_3P: -1.0,          # sn-Glycerol 3-phosphate
    Palmitate: -1.023,          # C16:0
    Stereate: -0.068,           # C18:0
    Water_c: -3.0,
    AMP: 3.0,
    Diphosphate: 3.0,
    Proton_c: 3.0,
    Phosphate: 1.0
})

reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

reaction = Reaction('TAG-Methanol')  
reaction.name = 'TAGs synthesis on Methanol'
reaction.subsystem = 'Lipid Metabolism'
reaction.lower_bound = 0.0
reaction.upper_bound =1000 # Irreversible

reaction.add_metabolites({
    ATP: -3.0,
    Margaroleate: -0.391,       # C17:1
    Heptadecanoate: -0.235,     # C17:0
    Tetradecanoate: -0.151,     # C14:0
    Pentadecanoate: -0.344,     # C15:0
    Pentadecenoic_acid: -0.032, # C15:1
    Palmitoleate: -0.399,       # C16:1
    Cis_Trans_Vaccenate: -0.472,# C18:1 (cis + trans)
    Glycerol_3P:-1.0,          # sn-Glycerol 3-phosphate
    Palmitate: -0.879,          # C16:0
    Stereate: -0.098,           # C18:0
    Water_c: -3.0,
    AMP: 3.0,
    Diphosphate: 3.0,
    Proton_c: 3.0,
    Phosphate: 1.0
})

reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

"""Transport Reactions"""

reaction = Reaction('Tr_acetate') 
reaction.name = 'Transport of Acetate symporter H'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.0
reaction.upper_bound = 1000.0  
reaction.add_metabolites({
    Acetate_e: -1.0,
    Proton_e: -1.0,
    Acetate_c: 1.0,
    Proton_c: 1.0,
})
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

reaction = Reaction('Tr_oxygen')  
reaction.name = 'Transport of Oxygen'
reaction.subsystem = 'Transport'
reaction.lower_bound = 0.0
reaction.upper_bound = 1000.0  
reaction.add_metabolites({
    Oxygen_c: 1.0,
    Oxygen_e: -1.0,
})
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

reaction = Reaction('Tr_CO2')  
reaction.name = 'Transport of CO2'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.0  
reaction.add_metabolites({
    CO2_c: 1.0,
    CO2_e: -1.0,
})
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

reaction = Reaction('Tr_Water')  
reaction.name = 'Transport of Water'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.0  
reaction.add_metabolites({
    Water_c: 1.0,
    Water_e: -1.0,
})
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

reaction = Reaction('Tr_H2')  
reaction.name = 'Transport of Hydrogen'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.0 
reaction.add_metabolites({
    H2_c: 1.0,
    H2_e: -1.0,
})
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

reaction = Reaction('Tr_Propanol') 
reaction.name = 'Transport of 1-Propanol, simple difussion'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.0  
reaction.add_metabolites({
    Propanol_c: 1.0,
    Propanol_e: -1.0,
})
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

reaction = Reaction('Tr_Butanoate') 
reaction.name = 'Transport of Butanoate symporter H'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.0  
reaction.add_metabolites({
    Butanoate_c: 1.0,
    Proton_c: 1.0,
    Proton_e:-1.0,
    Butanoate_e: -1.0,
})
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

reaction = Reaction('Tr_Valerate')  
reaction.name = 'Transport of Valerate symporter H'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.0 
reaction.add_metabolites({
    Valerate_c: 1.0,
    Proton_c: 1.0,
    Proton_e:-1.0,
    Valerate_e: -1.0,
})
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

reaction = Reaction('Tr_MeOH')  
reaction.name = 'Transport of Methanol simple diffusion'
reaction.subsystem = 'Transport'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.0  
reaction.add_metabolites({
    Methanol_c: 1.0,
    Methanol_e: -1.0,
})
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

reaction = Reaction('Tr_Polysaccharide')  
reaction.name = 'Transport of Polysaccharide '
reaction.subsystem = 'Unknown transporter'
reaction.lower_bound = -1000
reaction.upper_bound = 1000.0  
reaction.add_metabolites({
    Glycogen_e: 1.0,
    Glycogen_c: -1.0,
})
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

model.add_boundary(model.metabolites.get_by_id("Acetate_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("Water_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("Oxygen_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("Proton_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("CO2_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("Glucose_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("H2_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("Propanol_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("Butanoate_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("Valerate_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("Methanol_e") , type = "exchange")
model.add_boundary(model.metabolites.get_by_id("Glycogen_e") , type = "exchange")

model.reactions.get_by_id("EX_Acetate_e").lower_bound= 0
model.reactions.get_by_id("EX_Acetate_e").upper_bound= 0

model.reactions.get_by_id("EX_Glucose_e").lower_bound= -1.94
model.reactions.get_by_id("EX_Glucose_e").upper_bound=  0

model.reactions.get_by_id("EX_H2_e").lower_bound =0
model.reactions.get_by_id("EX_H2_e").upper_bound = 0

model.reactions.get_by_id("EX_Propanol_e").lower_bound = 0
model.reactions.get_by_id("EX_Propanol_e").upper_bound = 0

model.reactions.get_by_id("EX_Butanoate_e").lower_bound = 0
model.reactions.get_by_id("EX_Butanoate_e").upper_bound = 0

model.reactions.get_by_id("EX_Valerate_e").lower_bound = 0
model.reactions.get_by_id("EX_Valerate_e").upper_bound = 0

model.reactions.get_by_id("EX_Methanol_e").lower_bound = 0
model.reactions.get_by_id("EX_Methanol_e").upper_bound = 0

model.reactions.get_by_id("EX_Oxygen_e").lower_bound = -1000
model.reactions.get_by_id("EX_Oxygen_e").upper_bound = 1000

model.reactions.get_by_id("EX_CO2_e").lower_bound = -1000
model.reactions.get_by_id("EX_CO2_e").upper_bound = 1000

model.reactions.get_by_id("EX_Water_e").lower_bound = -1000
model.reactions.get_by_id("EX_Water_e").upper_bound = 1000

model.reactions.get_by_id("EX_Proton_e").lower_bound = -1000
model.reactions.get_by_id("EX_Proton_e").upper_bound = 1000


model.reactions.get_by_id("TAG-acetate").objective_coefficient= 0
model.reactions.get_by_id("TAG-Knallgas").objective_coefficient= 0
model.reactions.get_by_id("TAG-ButiricAcid").objective_coefficient= 0
model.reactions.get_by_id("TAG-1-Propanol").objective_coefficient= 0
model.reactions.get_by_id("TAG-Glucose").objective_coefficient= 1
model.reactions.get_by_id("TAG-ValericAcid").objective_coefficient= 0
model.reactions.get_by_id("TAG-Methanol").objective_coefficient= 0
model.reactions.get_by_id("Glycogen_synthase").objective_coefficient= 0
fba = model.optimize()

"""# Resultados"""

print(model.summary(fba))
print(fba.status)
