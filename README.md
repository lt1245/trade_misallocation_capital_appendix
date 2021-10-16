# trade_misallocation_capital_appendix
Online appendix for Trade, Misallocation, and Capital Market Integration

code_new.jl: Contains objects necessary for the steady states. All the packages must be installed, and at least 2 CPU-s must be available.

main_script.jl: Steady states evaluated. Calls code_new.jl

transition.jl: Contains objects necessary for transition paths. Assumes main_script.jl is evaluated, and case_final must be set appropriately, for values 1 to 8.

plot.jl: Assumes that transition.jl is evaluated for each case_final values
