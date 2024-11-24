from fsam_physics_model import calculate_equivalent_plastic_strain
import pandas as pd

def get_model_parameters(experiment='default'):
  """
  Gets the model parameters for the experiment, default is the example for 6061 in the original paper
  """
  johnson_cook_paper_data_6061 = { # from paper for 6061 aluminum 
    "Model": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
    "A (MPa)": [250, 293.4, 324.1, 250, 250, 250, 250, 335, 236.7, 275, 164],
    "B (MPa)": [79.7, 121.2, 113.8, 70, 70, 137, 209, 85, 41.2, 86, 211],
    "C": [0.0249, 0.002, 0.002, 0.001, 0.001, 0.0205, 0.001, 0.012, 0.0411, 0.0031, 0.0019],
    "n": [0.499, 0.23, 0.42, 0.499, 0.499, 0.499, 0.499, 0.1, 0.084, 0.39, 0.465],
    "m": [1.499, 1.34, 1.34, 1, 1.315, 1.499, 1.499, 1, 1.41, 1, 1.419],
    "Reference": [27, 28, 29, 30, 30, 30, 30, 31, 32, 33, 34],
  }
  jc_6061_data = pd.DataFrame(johnson_cook_paper_data_6061)
  numerical_columns = [col for col in jc_6061_data.columns if col not in ['Model', 'Reference']]
  jc_6061_data_means = jc_6061_data[numerical_columns].mean().to_dict()
    
  model_parameters = {}
  if experiment == 'default':
    feedstock_side_length = 12.7 # mm
    friction_force_radius = 9.525  # mm
    layer_thickness = 1.5  # mm
    
    model_parameters = {
      "A": jc_6061_data_means["A (MPa)"],
      "B": jc_6061_data_means["B (MPa)"],
      "C": jc_6061_data_means["C"],
      "n": jc_6061_data_means["n"],
      "m": jc_6061_data_means["m"],
      "friction_force_radius": friction_force_radius,
      "layer_thickness": layer_thickness,
      "feedstock_side_length": feedstock_side_length,
      "equivalent_plastic_strain":  calculate_equivalent_plastic_strain(feedstock_side_length, friction_force_radius, layer_thickness),
      "feed_velocity": 126,  # mm/min
      "temperature_initial": 25,  # deg C
      "temperature_feedstock_melting": 652,  # deg C
      "coulomb_friction_coefficient": 0.25,
      "effective_strain_rate": None,
      "reference_strain_rate": 1,  # According to the Johnson Cook model
      "temperature_deposition": 370,  # deg C
      "thermal_conductivity": 160  # W/m·K (from Figure 4)
    }

  elif experiment == 'mini':
    # Modified for smaller feed stock and tool size
    tool_diameter = 9.525 # mm
    feedstock_side_length = 3.175 # mm
    friction_force_radius  = (tool_diameter/2)/2 # mm. half of the tool radius
    layer_thickness = .5  # mm
    
    print(f'Tool diameter: {tool_diameter}mm')
    print(f'Feedstock side length: {feedstock_side_length}mm')
    print(f'Layer thickness: {layer_thickness}mm')
    print()
    model_parameters = {
      "tool_diameter": tool_diameter,
      "A": jc_6061_data_means["A (MPa)"],
      "B": jc_6061_data_means["B (MPa)"],
      "C": jc_6061_data_means["C"],
      "n": jc_6061_data_means["n"],
      "m": jc_6061_data_means["m"],
      "friction_force_radius": friction_force_radius,
      "layer_thickness": layer_thickness,
      "feedstock_side_length": feedstock_side_length,
      "equivalent_plastic_strain":  calculate_equivalent_plastic_strain(feedstock_side_length, friction_force_radius, layer_thickness),
      "feed_velocity": 126,  # mm/min
      "temperature_initial": 25,  # deg C
      "temperature_feedstock_melting": 652,  # deg C
      "coulomb_friction_coefficient": 0.25,
      "effective_strain_rate": None,
      "reference_strain_rate": 1,  # According to the Johnson Cook model
      "temperature_deposition": 370,  # deg C
      "thermal_conductivity": 160  # W/m·K (from Figure 4)
    }

  return model_parameters
