from fsam_physics_model import calculate_spindle_speed_given_deposition_temperature, calculate_normal_force
from parameters import get_model_parameters
import numpy as np

# Define the deposition temperatures that you want to evaluate the spindle speed at
temperature_deposition_array = [250, 300, 350, 400, 450] # for 6061 aluminum

# and their corresponding temeprature dependent coefficient
thermal_conductivity_array = [179, 172.5, 164, 155.5, 152.5] # for 6061 aluminum

model_parameters = get_model_parameters('mini') # modified parameters for smaller tool
print(f"Required normal force: {calculate_normal_force(np.pi*(model_parameters['tool_diameter']/2)**2)}N")
for i in range(len(temperature_deposition_array)):
  model_parameters['temperature_deposition'] = temperature_deposition_array[i]
  model_parameters['thermal_conductivity'] = thermal_conductivity_array[i]
  spindle_speed = calculate_spindle_speed_given_deposition_temperature(model_parameters)
  print(f'Spindle speed: {spindle_speed:.2f} rpm    Depostition temperature: {temperature_deposition_array[i]}°C ')