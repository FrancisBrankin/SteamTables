#Constant Volume Testing
#T_1, T_2, P_1, P_2, q_12, c_v, s_1, s_2
# if __name__ == "__main__":
#     inputs = 8
#     eqn_values = 5
#     for a in range(0,eqn_values):
#         values = [4,8,2,6,0.454822556]
#         #T_1, T_2, P_1, P_2, q_12, c_v, s_1, s_2
#         positions = [5,1,0,7,6]
#         values.pop(a)
#         positions.pop(a)
#         #print(str(values))
#         counter = 0
#         variables = [None] * inputs
#         for i in range(0,inputs):
#             if i in positions:
#                 variables[positions[counter]] = values[counter]
#                 #print(variables[positions[counter]])
#                 counter += 1

#             #else:
#                 #variables.append(None)

#         #print(variables)

#         test = Constant_Volume(variables[0],variables[1],variables[2],
#                 variables[3],variables[4],variables[5],variables[6],variables[7])
#         print(test.equation_finder()[0])

    #def __init__(self, T_1, T_2, P_1, P_2, q_12, c_v, s_1, s_2):

# #Constant Pressure Testing
# #     #T_1, T_2, v_1, v_2, q_12, w_12, c_v, c_p, R, s_1, s_2, p
# if __name__ == "__main__":
#     inputs = 12
#     eqn_values = 5
#     for a in range(0,eqn_values):
#         values = [4,8,2,6,0.4548225555]
#         #T_1, T_2, P_1, P_2, q_12, c_v, s_1, s_2
#         positions = [7,1,0,10,9]
#         values.pop(a)
#         positions.pop(a)
#         print(str(values))
#         counter = 0
#         variables = [None] * inputs
#         for i in range(0,inputs):
#             if i in positions:
#                 variables[positions[counter]] = values[counter]
#                 #print(variables[positions[counter]])
#                 counter += 1

#             else:
#                 variables.append(None)

#         #print(variables)


#         test = Constant_Pressure(variables[0],variables[1],variables[2],
#             variables[3],variables[4],variables[5],variables[6],variables[7],
#             variables[8],variables[9],variables[10],variables[11])
#         print(test.equation_finder()[0])

# #Constant Temperature Testing
# # T, v_1, v_2, P_1, P_2, q_12, w_12, c_v, c_p, R, s_1, s_2
# if __name__ == "__main__":
#     inputs = 12
#     eqn_values = 5
#     for a in range(0,eqn_values):
#         values = [41.58883083,5,6,8,2]
#         #T, v_1, v_2, P_1, P_2, q_12, w_12, c_v, c_p, R, s_1, s_2
#         positions = [5,0,9,2,1]
#         values.pop(a)
#         positions.pop(a)
#         print(str(values))
#         counter = 0
#         variables = [None] * inputs
#         for i in range(0,inputs):
#             if i in positions:
#                 variables[positions[counter]] = values[counter]
#                 #print(variables[positions[counter]])
#                 counter += 1

#             else:
#                 variables.append(None)

#         #print(variables)


#         test = Constant_Temperature(variables[0],variables[1],variables[2],
#             variables[3],variables[4],variables[5],variables[6],variables[7],
#             variables[8],variables[9],variables[10],variables[11])
#         print(test.equation_finder()[0])

# #Adiabatic Testing
# #T_1, T_2, v_1, v_2, P_1, P_2, w_12, c_v, c_p, R, s_1, s_2, u_1, u_2, gamma

# if __name__ == "__main__":
#     inputs = 15
#     eqn_values = 4
#     for a in range(0,eqn_values):
#         values = [-24,4,8,2]
#         #T_1, T_2, v_1, v_2, P_1, P_2, w_12, c_v, c_p, R, s_1, s_2, u_1, u_2, gamma
#         # 0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6  ,  7 ,  8 , 9, 10 , 11 , 12 , 13 ,  14
#         positions = [6,7,1,0]
#         values.pop(a)
#         positions.pop(a)
#         print(str(values))
#         counter = 0
#         variables = [None] * inputs
#         for i in range(0,inputs):
#             if i in positions:
#                 variables[positions[counter]] = values[counter]
#                 #print(variables[positions[counter]])
#                 counter += 1

#             else:
#                 variables.append(None)

#         #print(variables)


#         test = Adiabatic(variables[0],variables[1],variables[2],
#             variables[3],variables[4],variables[5],variables[6],variables[7],
#             variables[8],variables[9],variables[10],variables[11],variables[12],
#             variables[13],variables[14])
#         print(test.equation_finder()[0])

# #Polytropic Testing
# #P_1, P_2, v_1, v_2, T_1, T_2, s_1, s_2, w_12, q_12, c_v, R, n
# if __name__ == "__main__":
#     inputs = 13
#     eqn_values = 7
#     for a in range(0,eqn_values):
#         values = [1,0.5288678574,5,6,7,8,9]
#         #P_1, P_2, v_1, v_2, T_1, T_2, s_1, s_2, w_12, q_12, c_v, R, n
#         # 0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  8  ,  9  , 10 ,11, 12
#         positions = [7,6,10,11,12,4,5]
#         values.pop(a)
#         positions.pop(a)
#         print(str(values))
#         counter = 0
#         variables = [None] * inputs
#         for i in range(0,inputs):
#             if i in positions:
#                 variables[positions[counter]] = values[counter]
#                 #print(variables[positions[counter]])
#                 counter += 1

#             else:
#                 variables.append(None)

#         #print(variables)


#         test = Polytropic(variables[0],variables[1],variables[2],
#             variables[3],variables[4],variables[5],variables[6],variables[7],
#             variables[8],variables[9],variables[10],variables[11],variables[12])
#         print(test.equation_finder()[0])

# #Flow Processes Testing
# #q_12, w_12, h_2, h_1, C_1, C_2, z_1, z_2
# if __name__ == "__main__":
#     inputs = 8
#     eqn_values = 8
#     for a in range(0,eqn_values):
#         values = [15,0.69,1,2,3,4,5,6]
#         #q_12, w_12, h_1, h_2, C_1, C_2, z_1, z_2
#         positions = [0,1,2,3,4,5,6,7]
#         values.pop(a)
#         positions.pop(a)
#         print(str(values))
#         counter = 0
#         variables = [None] * inputs
#         for i in range(0,inputs):
#             if i in positions:
#                 variables[positions[counter]] = values[counter]
#                 #print(variables[positions[counter]])
#                 counter += 1

#             else:
#                 variables.append(None)

#         #print(variables)


#         test = Flow_Processes(variables[0],variables[1],variables[2],
#             variables[3],variables[4],variables[5],variables[6],variables[7])
#         print(test.equation_finder()[0])

# #Boilers Condesnors Heaters Coolers Testing
# #q_12, T_1, T_2, h_1, h_2, s_1, s_2, c_p
# if __name__ == "__main__":
#     inputs = 8
#     eqn_values = 4
#     for a in range(0,eqn_values):
#         values = [4,2,2,1]
#         #T_1, T_2, P_1, P_2, q_12, c_v, s_1, s_2
#         positions = [1,0,3,2]
#         values.pop(a)
#         positions.pop(a)
#         print(str(values))
#         counter = 0
#         variables = [None] * inputs
#         for i in range(0,inputs):
#             if i in positions:
#                 variables[positions[counter]] = values[counter]
#                 #print(variables[positions[counter]])
#                 counter += 1

#             else:
#                 variables.append(None)

#         #print(variables)


#         test = Boilers_Condensors_Heaters_Coolers(variables[0],variables[1],variables[2],
#             variables[3],variables[4],variables[5],variables[6],variables[7])
#         print(test.equation_finder()[0])

# #Nozzles Diffusers Testing
# #C_1, C_2, h_1, h_2, h_2i, T_1, T_2, T_2i, s_1, s_2, c_p, mu
# if __name__ == "__main__":
#     inputs = 12
#     eqn_values = 5
#     for a in range(0,eqn_values):
#         values = [3,4,5,1,0.33056934605]
#         #C_1, C_2, h_1, h_2, h_2i, T_1, T_2, T_2i, s_1, s_2, c_p, mu
#         # 0 ,  1 ,  2 ,  3 ,  4  ,  5 ,  6 ,  7  ,  8 ,  9 , 10 , 11
#         positions = [10,7,6,9,8]
#         values.pop(a)
#         positions.pop(a)
#         print(str(values))
#         counter = 0
#         variables = [None] * inputs
#         for i in range(0,inputs):
#             if i in positions:
#                 variables[positions[counter]] = values[counter]
#                 #print(variables[positions[counter]])
#                 counter += 1

#             else:
#                 variables.append(None)

#         #print(variables)


#         test = Nozzles_Diffusers(variables[0],variables[1],variables[2],
#             variables[3],variables[4],variables[5],variables[6],variables[7],
#             variables[8],variables[9],variables[10],variables[11])
#         print(test.equation_finder()[0])

# #Turbine Compressors Testing
# #w_12,h_1,h_2,h_2i,T_1,T_2,T_2i,s_1,s_2,c_p
# if __name__ == "__main__":
#     inputs = 10
#     eqn_values = 4
#     for a in range(0,eqn_values):
#         values = [4,2,2,1]
#         #T_1, T_2, P_1, P_2, q_12, c_v, s_1, s_2
#         positions = [1,0,3,2]
#         values.pop(a)
#         positions.pop(a)
#         print(str(values))
#         counter = 0
#         variables = [None] * inputs
#         for i in range(0,inputs):
#             if i in positions:
#                 variables[positions[counter]] = values[counter]
#                 #print(variables[positions[counter]])
#                 counter += 1

#             else:
#                 variables.append(None)

#         #print(variables)


#     test = Turbine_Compressors(variables[0],variables[1],variables[2],
#             variables[3],variables[4],variables[5],variables[6],variables[7],
#             variables[8],variables[9])
#     print(test.equation_finder()[0])

# #Throttles Testing
# #h_1,h_2
# if __name__ == "__main__":
#     inputs = 2
#     eqn_values = 4
#     for a in range(0,eqn_values):
#         values = [4,2,2,1]
#         #T_1, T_2, P_1, P_2, q_12, c_v, s_1, s_2
#         positions = [1,0,3,2]
#         values.pop(a)
#         positions.pop(a)
#         print(str(values))
#         counter = 0
#         variables = [None] * inputs
#         for i in range(0,inputs):
#             if i in positions:
#                 variables[positions[counter]] = values[counter]
#                 #print(variables[positions[counter]])
#                 counter += 1

#             else:
#                 variables.append(None)

#         #print(variables)


#     test = Throttles(variables[0],variables[1])
#     print(test.equation_finder()[0])

