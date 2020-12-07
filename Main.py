from MES import grid, global_data

grid.print_grid_data()
# grid.print_global_matrix()
grid.plot_grid()
# grid.to_file()
tmp = 0
tmp2 = []
for i in range(global_data.nE):
    tmp = grid.ELEM[i].boundary_condition()
    tmp2.append(tmp[0])
print(grid.P_global(tmp2))
