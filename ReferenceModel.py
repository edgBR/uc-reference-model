########################################################################################################
# a basic (thermal) unit commitment model, drawn from:                                                 #
# A Computationally Efficient Mixed-Integer Linear Formulation for the Thermal Unit Commitment Problem #
# Miguel Carrion and Jose M. Arroyo                                                                    #
# IEEE Transactions on Power Systems, Volume 21, Number 3, August 2006.                                #
########################################################################################################

from coopr.pyomo import *

model = Model()

#
# Parameters
#

##########################################################
# string indentifiers for the set of thermal generators. #
##########################################################

model.ThermalGenerators = Set()

###################################################
# the number of time periods under consideration, #
# in addition to the corresponding set.           #
###################################################

model.NumTimePeriods = Param(within=PositiveIntegers)

model.TimePeriods = RangeSet(1, model.NumTimePeriods)

#################################################################
# the global system demand, for each time period. units are MW. #
#################################################################

model.Demand = Param(model.TimePeriods, within=NonNegativeReals)

##################################################################
# the global system reserve, for each time period. units are MW. #
##################################################################

model.ReserveRequirement = Param(model.TimePeriods, within=NonNegativeReals, default=0.0)

####################################################################################
# minimum and maximum generation levels, for each thermal generator. units are MW. #
# could easily be specified on a per-time period basis, but are not currently.     #
####################################################################################

model.MinimumPowerOutput = Param(model.ThermalGenerators, within=NonNegativeReals, default=0.0)

def maximum_power_output_validator(m, v, g):
   return v >= value(m.MinimumPowerOutput[g])

model.MaximumPowerOutput = Param(model.ThermalGenerators, within=NonNegativeReals, validate=maximum_power_output_validator)

#################################################
# generator ramp up/down rates. units are MW/h. #
#################################################

# limits for normal time periods
model.NominalRampUpLimit = Param(model.ThermalGenerators, within=NonNegativeReals)
model.NominalRampDownLimit = Param(model.ThermalGenerators, within=NonNegativeReals)

# limits for time periods in which generators are brought on or off-line. 
# must be no less than the generator minimum output. 
def at_least_generator_minimum_output_validator(m, v, g):
   return v >= m.MinimumPowerOutput[g]

model.StartupRampLimit = Param(model.ThermalGenerators, within=NonNegativeReals, validate=at_least_generator_minimum_output_validator)
model.ShutdownRampLimit = Param(model.ThermalGenerators, within=NonNegativeReals, validate=at_least_generator_minimum_output_validator)

##########################################################################################################
# the minimum number of time periods that a generator must be on-line (off-line) once brought up (down). #
##########################################################################################################

model.MinimumUpTime = Param(model.ThermalGenerators, within=NonNegativeIntegers)
model.MinimumDownTime = Param(model.ThermalGenerators, within=NonNegativeIntegers)

#############################################
# unit on state at t=0 (initial condition). #
#############################################

# if positive, the number of hours prior to (and including) t=0 that the unit has been on.
# if negative, the number of hours prior to (and including) t=0 that the unit has been off.
# the value cannot be 0, by definition.

def t0_state_nonzero_validator(m, v, g):
    return v != 0

model.UnitOnT0State = Param(model.ThermalGenerators, within=Integers, validate=t0_state_nonzero_validator)

def t0_unit_on_rule(m, g):
    return value(m.UnitOnT0State[g]) >= 1

model.UnitOnT0 = Param(model.ThermalGenerators, within=Binary, initialize=t0_unit_on_rule)

#######################################################################################
# the number of time periods that a generator must initally on-line (off-line) due to #
# its minimum up time (down time) constraint.                                         #
#######################################################################################

def initial_time_periods_online_rule(m, g):
   if not value(m.UnitOnT0[g]):
      return 0
   else:
      return min(value(m.NumTimePeriods), \
                 max(0, \
                     value(m.MinimumUpTime[g]) - value(m.UnitOnT0State[g])))

model.InitialTimePeriodsOnLine = Param(model.ThermalGenerators, within=NonNegativeIntegers, initialize=initial_time_periods_online_rule)

def initial_time_periods_offline_rule(m, g):
   if value(m.UnitOnT0[g]):
      return 0
   else:
      return min(value(m.NumTimePeriods), \
                 max(0, \
                     value(m.MinimumDownTime[g]) + value(m.UnitOnT0State[g]))) # m.UnitOnT0State is negative if unit is off

model.InitialTimePeriodsOffLine = Param(model.ThermalGenerators, within=NonNegativeIntegers, initialize=initial_time_periods_offline_rule)

####################################################################
# generator power output at t=0 (initial condition). units are MW. #
####################################################################

model.PowerGeneratedT0 = Param(model.ThermalGenerators, within=NonNegativeReals)

##################################################################################################################
# production cost coefficients (for the quadratic) a0=constant, a1=linear coefficient, a2=quadratic coefficient. #
##################################################################################################################

model.ProductionCostA0 = Param(model.ThermalGenerators, within=NonNegativeReals) # units are $/hr (or whatever the time unit is).
model.ProductionCostA1 = Param(model.ThermalGenerators, within=NonNegativeReals) # units are $/MWhr.
model.ProductionCostA2 = Param(model.ThermalGenerators, within=NonNegativeReals) # units are $/(MWhr^2).

##############################################################################################
# number of pieces in the linearization of each generator's quadratic cost production curve. #
##############################################################################################

model.NumGeneratorCostCurvePieces = Param(within=PositiveIntegers, default=2)

#######################################################################
# points for piecewise linearization of power generation cost curves. #
#######################################################################

# maps a (generator, time-index) pair to a list of points defining the piecewise cost linearization breakpoints.
# the time index is redundant, but required. in the Piecewise construct, the breakpoints must be indexed the 
# same as the Piecewise construct itself.

model.PowerGenerationPiecewisePoints = {} 

def power_generation_piecewise_points_rule(model, g, t):
    min_power = value(model.MinimumPowerOutput[g])
    max_power = value(model.MaximumPowerOutput[g])
    n = value(model.NumGeneratorCostCurvePieces)
    width = (max_power - min_power) / float(n)
    model.PowerGenerationPiecewisePoints[g, t] = [min_power + i*width for i in xrange(0,n+1)]

model.CreatePowerGenerationPiecewisePoints = BuildAction(model.ThermalGenerators * model.TimePeriods, rule=power_generation_piecewise_points_rule)

# a function for use in piecewise linearization of the cost function. 
def production_cost_function(model, g, t, x):
    return value(model.ProductionCostA0[g]) + value(model.ProductionCostA1[g])*x + value(model.ProductionCostA2[g])*x*x

##################################################################################
# shutdown cost for each generator. in the literature, these are often set to 0. #
##################################################################################

model.ShutdownCostCoefficient = Param(model.ThermalGenerators, within=NonNegativeReals, default=0.0) # units are $.

#
# Variables 
#

######################
# decision variables #
######################

# indicator variables for each generator, at each time period.
model.UnitOn = Var(model.ThermalGenerators, model.TimePeriods, within=Binary)

# amount of power produced by each generator, at each time period.
model.PowerGenerated = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals)

# maximum power output for each generator, at each time period.
model.MaximumPowerAvailable = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals)

###################
# cost components #
###################

# production cost associated with each generator, for each time period.
model.ProductionCost = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals)

# startup and shutdown costs for each generator, each time period.
model.StartupCost = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals)
model.ShutdownCost = Var(model.ThermalGenerators, model.TimePeriods, within=NonNegativeReals)

# cost over all generators, for all time periods.
model.TotalProductionCost = Var(within=NonNegativeReals)

# all other overhead / fixed costs, e.g., associated with startup and shutdown.
model.TotalFixedCost = Var(within=NonNegativeReals)

#
# Constraints
#

# meet the demand at each time period.
# encodes Constraint 2 in Carrion and Arroyo.
def production_equals_demand_rule(m, t):
   return sum(m.PowerGenerated[g, t] for g in m.ThermalGenerators) == m.Demand[t]

model.ProductionEqualsDemand = Constraint(model.TimePeriods, rule=production_equals_demand_rule)

# ensure there is sufficient maximal power output available to meet both the 
# demand and the spinning reserve requirements in each time period.
# encodes Constraint 3 in Carrion and Arroyo.
def enforce_reserve_requirements_rule(m, t):
   return sum(m.MaximumPowerAvailable[g, t] for g in m.ThermalGenerators) >= m.Demand[t] + m.ReserveRequirement[t]

model.EnforceReserveRequirements = Constraint(model.TimePeriods, rule=enforce_reserve_requirements_rule)

############################################
# generation limit and ramping constraints #
############################################

# enforce the generator power output limits on a per-period basis.
# the maximum power available at any given time period is dynamic,
# bounded from above by the maximum generator output.

# the following three constraints encode Constraints 16 and 17 defined in Carrion and Arroyo.

# NOTE: The expression below is what we really want - however, due to a pyomo bug, we have to split it into two constraints:
# m.MinimumPowerOutput[g] * m.UnitOn[g, t] <= m.PowerGenerated[g,t] <= m.MaximumPowerAvailable[g, t]
# When fixed, merge back parts "a" and "b", leaving two constraints.

def enforce_generator_output_limits_rule_part_a(m, g, t):
   return m.MinimumPowerOutput[g] * m.UnitOn[g, t] <= m.PowerGenerated[g,t]

model.EnforceGeneratorOutputLimitsPartA = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_generator_output_limits_rule_part_a)

def enforce_generator_output_limits_rule_part_b(m, g, t):
   return m.PowerGenerated[g,t] <= m.MaximumPowerAvailable[g, t]

model.EnforceGeneratorOutputLimitsPartB = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_generator_output_limits_rule_part_b)

def enforce_generator_output_limits_rule_part_c(m, g, t):
   return m.MaximumPowerAvailable[g,t] <= m.MaximumPowerOutput[g] * m.UnitOn[g, t]

model.EnforceGeneratorOutputLimitsPartC = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_generator_output_limits_rule_part_c)

# impose upper bounds on the maximum power available for each generator in each time period, 
# based on standard and start-up ramp limits.

# the following constraint encodes Constraint 18 defined in Carrion and Arroyo.

def enforce_max_available_ramp_up_rates_rule(m, g, t):
   # 4 cases, split by (t-1, t) unit status (RHS is defined as the delta from m.PowerGenerated[g, t-1])
   # (0, 0) - unit staying off:   RHS = maximum generator output (degenerate upper bound due to unit being off) 
   # (0, 1) - unit switching on:  RHS = startup ramp limit 
   # (1, 0) - unit switching off: RHS = standard ramp limit minus startup ramp limit plus maximum power generated (degenerate upper bound due to unit off)
   # (1, 1) - unit staying on:    RHS = standard ramp limit
   if t == 1:
      return m.MaximumPowerAvailable[g, t] <= m.PowerGeneratedT0[g] + \
                                              m.NominalRampUpLimit[g] * m.UnitOnT0[g] + \
                                              m.StartupRampLimit[g] * (m.UnitOn[g, t] - m.UnitOnT0[g]) + \
                                              m.MaximumPowerOutput[g] * (1 - m.UnitOn[g, t])
   else:
      return m.MaximumPowerAvailable[g, t] <= m.PowerGenerated[g, t-1] + \
                                              m.NominalRampUpLimit[g] * m.UnitOn[g, t-1] + \
                                              m.StartupRampLimit[g] * (m.UnitOn[g, t] - m.UnitOn[g, t-1]) + \
                                              m.MaximumPowerOutput[g] * (1 - m.UnitOn[g, t])

model.EnforceMaxAvailableRampUpRates = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_max_available_ramp_up_rates_rule)

# the following constraint encodes Constraint 19 defined in Carrion and Arroyo.

def enforce_max_available_ramp_down_rates_rule(m, g, t):
   # 4 cases, split by (t, t+1) unit status
   # (0, 0) - unit staying off:   RHS = 0 (degenerate upper bound)
   # (0, 1) - unit switching on:  RHS = maximum generator output minus shutdown ramp limit (degenerate upper bound) - this is the strangest case.
   # (1, 0) - unit switching off: RHS = shutdown ramp limit
   # (1, 1) - unit staying on:    RHS = maximum generator output (degenerate upper bound)
   if t == value(m.NumTimePeriods):
      return Constraint.Skip
   else:
      return m.MaximumPowerAvailable[g, t] <= \
             m.MaximumPowerOutput[g] * m.UnitOn[g, t+1] + \
             m.ShutdownRampLimit[g] * (m.UnitOn[g, t] - m.UnitOn[g, t+1])

model.EnforceMaxAvailableRampDownRates = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_max_available_ramp_down_rates_rule)

# the following constraint encodes Constraint 20 defined in Carrion and Arroyo.

def enforce_ramp_down_limits_rule(m, g, t):
   # 4 cases, split by (t-1, t) unit status: 
   # (0, 0) - unit staying off:   RHS = maximum generator output (degenerate upper bound)
   # (0, 1) - unit switching on:  RHS = standard ramp-down limit minus shutdown ramp limit plus maximum generator output - this is the strangest case.
   # (1, 0) - unit switching off: RHS = shutdown ramp limit 
   # (1, 1) - unit staying on:    RHS = standard ramp-down limit 
   if t == 1:
      return m.PowerGeneratedT0[g] - m.PowerGenerated[g, t] <= \
             m.NominalRampDownLimit[g] * m.UnitOn[g, t] + \
             m.ShutdownRampLimit[g] * (m.UnitOnT0[g] - m.UnitOn[g, t]) + \
             m.MaximumPowerOutput[g] * (1 - m.UnitOnT0[g])                
   else:
      return m.PowerGenerated[g, t-1] - m.PowerGenerated[g, t] <= \
             m.NominalRampDownLimit[g] * m.UnitOn[g, t] + \
             m.ShutdownRampLimit[g] * (m.UnitOn[g, t-1] - m.UnitOn[g, t]) + \
             m.MaximumPowerOutput[g] * (1 - m.UnitOn[g, t-1])             

model.EnforceNominalRampDownLimits = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_ramp_down_limits_rule)

#############################################
# constraints for computing cost components #
#############################################

# compute the per-generator, per-time period production costs. this is a "simple" piecewise linear construct.
# the first argument to piecewise is the index set. the second and third arguments are respectively the input and output variables. 
model.ComputeProductionCosts = Piecewise(model.ThermalGenerators * model.TimePeriods, model.ProductionCost, model.PowerGenerated, pw_pts=model.PowerGenerationPiecewisePoints, f_rule=production_cost_function, pw_constr_type='LB')

# compute the total production costs, across all generators and time periods.
def compute_total_production_cost_rule(m):
   return m.TotalProductionCost == sum(m.ProductionCost[g, t] for g in m.ThermalGenerators for t in m.TimePeriods)

model.ComputeTotalProductionCost = Constraint(rule=compute_total_production_cost_rule)

# compute the per-generator, per-time period shutdown costs.
def compute_shutdown_costs_rule(m, g, t):
   if t is 1:
      return m.ShutdownCost[g, t] >= m.ShutdownCostCoefficient[g] * (m.UnitOnT0[g] - m.UnitOn[g, t])
   else:
      return m.ShutdownCost[g, t] >= m.ShutdownCostCoefficient[g] * (m.UnitOn[g, t-1] - m.UnitOn[g, t])

model.ComputeShutdownCosts = Constraint(model.ThermalGenerators, model.TimePeriods, rule=compute_shutdown_costs_rule)

# compute the total startup and shutdown costs, across all generators and time periods.
def compute_total_fixed_cost_rule(m):
   return m.TotalFixedCost == sum(m.StartupCost[g, t] + m.ShutdownCost[g, t] for g in m.ThermalGenerators for t in m.TimePeriods)

model.ComputeTotalFixedCost = Constraint(rule=compute_total_fixed_cost_rule)

#######################
# up-time constraints #
#######################

# constraint due to initial conditions.
def enforce_up_time_constraints_initial(m, g):
   if value(m.InitialTimePeriodsOnLine[g]) is 0:
      return Constraint.Skip
   return sum((1 - m.UnitOn[g, t]) for g in m.ThermalGenerators for t in m.TimePeriods if t <= value(m.InitialTimePeriodsOnLine[g])) == 0.0

model.EnforceUpTimeConstraintsInitial = Constraint(model.ThermalGenerators, rule=enforce_up_time_constraints_initial)

# constraint for each time period after that not involving the initial condition.
def enforce_up_time_constraints_subsequent(m, g, t):
   if t <= value(m.InitialTimePeriodsOnLine[g]):
      # handled by the EnforceUpTimeConstraintInitial constraint.
      return Constraint.Skip
   elif t <= (value(m.NumTimePeriods) - value(m.MinimumUpTime[g]) + 1):
      # the right-hand side terms below are only positive if the unit was off in the previous time period but on in this one =>
      # the value is the minimum number of subsequent consecutive time periods that the unit is required to be on.
      if t is 1:
         return sum(m.UnitOn[g, n] for n in m.TimePeriods if n >= t and n <= (t + value(m.MinimumUpTime[g]) - 1)) >= \
                (m.MinimumUpTime[g] * (m.UnitOn[g, t] - m.UnitOnT0[g]))
      else:
         return sum(m.UnitOn[g, n] for n in m.TimePeriods if n >= t and n <= (t + value(m.MinimumUpTime[g]) - 1)) >= \
                (m.MinimumUpTime[g] * (m.UnitOn[g, t] - m.UnitOn[g, t-1]))
   else:
      # handle the final (MinimumUpTime[g] - 1) time periods - if a unit is started up in 
      # this interval, it must remain on-line until the end of the time span.
      if t == 1: # can happen when small time horizons are specified
         return sum((m.UnitOn[g, n] - (m.UnitOn[g, t] - m.UnitOnT0[g])) for n in m.TimePeriods if n >= t) >= 0.0
      else:
         return sum((m.UnitOn[g, n] - (m.UnitOn[g, t] - m.UnitOn[g, t-1])) for n in m.TimePeriods if n >= t) >= 0.0

model.EnforceUpTimeConstraintsSubsequent = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_up_time_constraints_subsequent)

#########################
# down-time constraints #
#########################

# constraint due to initial conditions.
def enforce_down_time_constraints_initial(m, g):
   if value(m.InitialTimePeriodsOffLine[g]) is 0: 
      return Constraint.Skip
   return sum(m.UnitOn[g, t] for g in m.ThermalGenerators for t in m.TimePeriods if t <= value(m.InitialTimePeriodsOffLine[g])) == 0.0

model.EnforceDownTimeConstraintsInitial = Constraint(model.ThermalGenerators, rule=enforce_down_time_constraints_initial)

# constraint for each time period after that not involving the initial condition.
def enforce_down_time_constraints_subsequent(m, g, t):
   if t <= value(m.InitialTimePeriodsOffLine[g]):
      # handled by the EnforceDownTimeConstraintInitial constraint.
      return Constraint.Skip
   elif t <= (value(m.NumTimePeriods) - value(m.MinimumDownTime[g]) + 1):
      # the right-hand side terms below are only positive if the unit was off in the previous time period but on in this one =>
      # the value is the minimum number of subsequent consecutive time periods that the unit is required to be on.
      if t is 1:
         return sum((1 - m.UnitOn[g, n]) for n in m.TimePeriods if n >= t and n <= (t + value(m.MinimumDownTime[g]) - 1)) >= \
                (m.MinimumDownTime[g] * (m.UnitOnT0[g] - m.UnitOn[g, t]))
      else:
         return sum((1 - m.UnitOn[g, n] for n in m.TimePeriods if n >= t and n <= (t + value(m.MinimumDownTime[g]) - 1))) >= \
                (m.MinimumDownTime[g] * (m.UnitOn[g, t-1] - m.UnitOn[g, t]))
   else:
      # handle the final (MinimumDownTime[g] - 1) time periods - if a unit is shut down in
      # this interval, it must remain off-line until the end of the time span.
      if t == 1: # can happen when small time horizons are specified
         return sum(((1 - m.UnitOn[g, n]) - (m.UnitOnT0[g] - m.UnitOn[g, t])) for n in m.TimePeriods if n >= t) >= 0.0
      else:
         return sum(((1 - m.UnitOn[g, n]) - (m.UnitOn[g, t-1] - m.UnitOn[g, t])) for n in m.TimePeriods if n >= t) >= 0.0

model.EnforceDownTimeConstraintsSubsequent = Constraint(model.ThermalGenerators, model.TimePeriods, rule=enforce_down_time_constraints_subsequent)

#
# Objectives
#

def total_cost_objective_rule(m):
   return m.TotalProductionCost + m.TotalFixedCost

model.TotalCostObjective = Objective(rule=total_cost_objective_rule, sense=minimize)
