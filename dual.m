function hat_psi = dual(hat_phi, H)

hat_psi = H*hat_phi;
hat_psi = hat_psi./vecnorm(hat_psi);