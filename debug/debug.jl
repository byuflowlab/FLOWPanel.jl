A_mu_mu = view(solver.G, 1:body.ncells, 1:body.ncells)
A_w_mu = view(solver.G, body.ncells+1:body.ncells+body.nsheddings, 1:body.ncells)
A_w_w = view(solver.G, body.ncells+1:body.ncells+body.nsheddings, body.ncells+1:body.ncells+body.nsheddings)
A_mu_w = view(solver.G, 1:body.ncells, body.ncells+1:body.ncells+body.nsheddings)
