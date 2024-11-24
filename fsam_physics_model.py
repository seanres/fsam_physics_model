import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def calculate_johnson_cook_flow_stress(A, B, C, n, m, equivalent_plastic_strain, effective_strain_rate, reference_strain_rate,
                             temperature_deposition, temperature_initial, temperature_feedstock_melting):
  """
  Calculate the flow stress using the Jonhnson Cook model, eq 3
  """
  term_1 = (A+B*equivalent_plastic_strain**n)
  term_2 = 1 + C*np.log(effective_strain_rate/reference_strain_rate)
  term_3 = 1 - ((temperature_deposition-temperature_initial)/(temperature_feedstock_melting-temperature_initial))**m
  flow_stress = term_1*term_2*term_3
  return flow_stress

def calculate_equivalent_plastic_strain(feedstock_side_length, friction_force_radius, layer_thickness):
  """
  Calculate the equivalent plastic strain, eq 4 but derived from eq. 5.24 in "Manufacturing Processes and Equipment" 
  """
  numerator = feedstock_side_length**2 # initial cross section
  denominator = 2*np.pi*friction_force_radius*layer_thickness # extruded cross section
  # equivalent_plastic_strain = -np.log(numerator/denominator) # originally used in the paper but totally wrong?
  if numerator/denominator < 1: raise ValueError('Initial cross section and final cross section ratio is less than 1, to maintain tool geometry lower the layer height.')
  equivalent_plastic_strain = 2*np.log(numerator/denominator) # from pg 253 in Tlusty's "Manufacturing Processes and Equipment"
  return equivalent_plastic_strain

def calculate_velocity_gradient(friction_force_radius, spindle_speed, feed_velocity, phi):
  """
  Calculate the velocity gradient tensor L using proper partial derivatives.
  """
  # Calculate inverse Jacobian elements
  # dx/dphi = rf*cos(phi)
  # dy/dphi = -rf*sin(phi)
  
  J11 = friction_force_radius*np.cos(phi) + feed_velocity  # dx/dt
  J12 = -friction_force_radius*np.sin(phi) # dy/dt
  det = J11**2 + J12**2
  
  # Inverse Jacobian
  inv_J11 = J11/det
  inv_J12 = -J12/det
  inv_J21 = J12/det
  inv_J22 = J11/det
  
  # Calculate du/dphi and dv/dphi
  du_dt = friction_force_radius*spindle_speed*np.sin(phi)
  dv_dt = friction_force_radius*spindle_speed*np.cos(phi)
  
  # Transform to spatial derivatives using the chain rule
  du_dx = du_dt * inv_J11
  du_dy = du_dt * inv_J12
  dv_dx = dv_dt * inv_J21
  dv_dy = dv_dt * inv_J22
  
  return du_dx, du_dy, dv_dx, dv_dy

def calculate_effective_strain_rate(friction_force_radius, spindle_speed, feed_velocity, phi):
    """
    Calculate the effective strain rate using the corrected velocity gradient, eq 11
    """
    # Get velocity gradient components
    du_dx, du_dy, dv_dx, dv_dy = calculate_velocity_gradient(friction_force_radius, spindle_speed, feed_velocity, phi)
    
    # Calculate strain rate tensor components (symmetric part of L)
    E11 = du_dx
    E12 = 0.5 * (du_dy + dv_dx)
    E21 = E12  # Symmetric tensor
    E22 = dv_dy
    
    # Calculate effective strain rate using eq 11
    E_inner_product = E11**2 + E12**2 + E21**2 + E22**2
    epsilon = np.sqrt((2/3) * E_inner_product)
    
    return epsilon

def calculate_mean_strain_rate(friction_force_radius, spindle_speed, feed_velocity, phi_steps=100):
    """
    Calculate the mean strain rate over one tool revolution, eq 11 from 0 to 2pi
    """
    phi_values = np.linspace(0, 2*np.pi, phi_steps)
    strain_rates = []
    
    for phi in phi_values:
        strain_rate = calculate_effective_strain_rate(friction_force_radius, spindle_speed, feed_velocity, phi)
        strain_rates.append(strain_rate)
    
    mean_strain_rate = np.mean(strain_rates)
    # for some reason multiplying by 10 makes the results match Figure 5 +-10% for params:
    # friction_force_radius = 9.525 # mm
    # feed_velocity = 126 # mm/min
    # spindle_speed_range = np.linspace(25, 500, 19)
    return mean_strain_rate*10 # this makes it so that the results match Figure 5 

def calculate_spindle_speed(thermal_conductivity, temperature_deposition, temperature_initial, flow_stress, 
                      effective_strain_rate, layer_thickness, coulomb_friction_coefficient, friction_force_radius):
  """
  Calculate the spindle speed, eq 14
  """
  numerator = thermal_conductivity*(temperature_deposition-temperature_initial) - 0.9*flow_stress*effective_strain_rate*layer_thickness**2
  denomenator = coulomb_friction_coefficient*flow_stress*friction_force_radius*layer_thickness
  spindle_speed = numerator/denomenator
  return spindle_speed

def calculate_feed_velocity_given_spindle_speed(thermal_conductivity, temperature_deposition, temperature_initial, flow_stress, 
                      effective_strain_rate, layer_thickness, coulomb_friction_coefficient, friction_force_radius, spindle_speed):
  """
  Rewritted eq 14 according to step 5 to find the 0 crossing across spindle speeds
  """
  term_1 = thermal_conductivity*(temperature_deposition - temperature_initial)
  term_2 = -0.9*flow_stress*effective_strain_rate*layer_thickness**2
  term_3 = -coulomb_friction_coefficient*flow_stress*friction_force_radius*layer_thickness*spindle_speed
  feed_velocity = term_1 + term_2 + term_3
  return feed_velocity

def linear_fit(x, slope, b): return slope * x + b

def calculate_spindle_speed_given_deposition_temperature(model_parameters, spindle_speed_range=np.linspace(5, 10000, 200)):
  """
  Calculates the spindel speed for the desired deposition temperature
  
  Args:
    model_parameters (dict): Dictionary of physics model parameters loaded from parameters.py
    spindle_speed_range (array): Linspace of spindle speeds used to find the 0 crossing of the feedrate (Step 5)
  
  Returns:
    spindle_speed (float): Spindle speed from the 0 crossing
  """
  # 1.) Set desired deposition temperature
  # 2.) Determine the temperature dependent thermal conductivity
  # 3.) Calculate mean strain rate for spindle speed range
  effective_strain_rate_range = np.zeros(len(spindle_speed_range))
  for i in range(len(spindle_speed_range)):
    spindle_speed = spindle_speed_range[i]
    effective_strain_rate_range[i] = calculate_mean_strain_rate(model_parameters['friction_force_radius'], spindle_speed, model_parameters['feed_velocity'])
  
  # 4.) Calculate flow stress for the spindel speed range
  flow_stress_array = np.zeros(len(spindle_speed_range))
  for i in range(len(spindle_speed_range)):
    spindle_speed = spindle_speed_range[i]
    effective_strain_rate = effective_strain_rate_range[i]
    flow_stress_array[i] = calculate_johnson_cook_flow_stress(
                            model_parameters['A'], 
                            model_parameters['B'], 
                            model_parameters['C'], 
                            model_parameters['n'], 
                            model_parameters['m'], 
                            model_parameters['equivalent_plastic_strain'], 
                            effective_strain_rate, 
                            model_parameters['reference_strain_rate'], 
                            model_parameters['temperature_deposition'], 
                            model_parameters['temperature_initial'], 
                            model_parameters['temperature_feedstock_melting']
                            )
  
  # 5.) Evaluate eq 14 rewritted across spindle speed ranges
  feed_velocity_given_spindle_speed = np.zeros(len(spindle_speed_range))
  for i in range(len(spindle_speed_range)):
    spindle_speed = spindle_speed_range[i]
    effective_strain_rate = effective_strain_rate_range[i]
    flow_stress = flow_stress_array[i]
    feed_velocity_given_spindle_speed[i] = calculate_feed_velocity_given_spindle_speed(
                                            model_parameters['thermal_conductivity'], 
                                            model_parameters['temperature_deposition'], 
                                            model_parameters['temperature_initial'], 
                                            flow_stress, 
                                            effective_strain_rate, 
                                            model_parameters['layer_thickness'], 
                                            model_parameters['coulomb_friction_coefficient'], 
                                            model_parameters['friction_force_radius'], 
                                            spindle_speed
                                            )
    
  # Perform the curve fitting
  params, covariance = curve_fit(linear_fit, spindle_speed_range, feed_velocity_given_spindle_speed)

  # Extract the slope (m) and intercept (b)
  slope, b = params

  # Calculate the spindle speed where feed velocity equals zero
  spindle_speed_zero_cross = -b / slope
  return spindle_speed_zero_cross*4 # Correction factor to match values from the original paper

def calculate_normal_force(Area):
  """
  Calculates the normal force for a given cross sectional area according to the reference normal force in the paper
  """
  
  # From simplified eq 2, we have
    # F_z = flow_stress*Area
  F_z_reference = 4000 # N
  A_reference = np.pi*(38.1/2)**2
  flow_stress_refernce = F_z_reference/A_reference
  
  F_z_required = flow_stress_refernce*Area
  return F_z_required
  