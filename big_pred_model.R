model{
  ###########################################
  # Data that needs to be supplied to model #
  ###########################################
  #
  # 1. pred_psi_design_matrix: nsite by n_parm_pred, first column should be 1's for intercept.
  # 2. prey_psi_design_matrix: nsite by n_parm_prey, first column should be 1's for intercept,
  #      currently assumes that every covariate you put in here has an interaction with
  #      predator presence. Could be modified to change that pretty easily.
  # 3. pred_y: nsite by npred by nseason array. Each element is the number of secondary sample units the
  #              predator was observed at a site during a seasonal survey. If a survey
  #              did not happen at that site for a given season it should be NA.
  # 4. prey_y: nsite by nprey by nseason array. Each element is the number of secondary sample units the
  #              prey was observed at a site during a seasonal survey. If a survey
  #              did not happen at that site for a given season it should be NA.
  # 5. pred_rho_design_matrix: nsite by n_parm_pred, first column should be 1's for intercept.
  # 6. prey_rho_design_matrix: nsite by n_parm_prey, first column should be 1's for intercept.
  # 7. J: nsite by nseason matrix. Each element is the number of secondary sample units during a primary
  #         sampling period at a given site. If no sampling occured, that cell element is 0.
  #
  #################################################
  # Indices that need to be supplied to the model #
  #################################################
  # 1. npred: scalar. The number of predator.
  # 2. nsite: scalar. The number of unique sites across the study.
  # 3. nseason: The number of primary sampling periods.
  # 4. nparm_prey_psi: The number of prey parameters for psi (i.e., columns in the prey_psi_design_matrix).
  # 5. nparm_prey_rho: The number of prey parameters for rho (i.e., columns in the prey_rho_design_matrix).
  # 6. nprey: The number of prey species.
  # 7. nparm_pred_psi: The number of predator parameters for psi (i.e., columns in the pred_psi_design_matrix).
  # 8. nparm_pred_rho: The number of predator parameters for rho (i.e., columns in the pred_rho_design_matrix).
  
# Latent state dominant species (ds) model
  for(ds in 1:npred){
    for(site in 1:nsite){
      # logit linear predictor for predators, season 1
      logit(pred_psi[site,ds,1]) <- inprod(
        pred_beta[ds,],
        pred_psi_design_matrix[site,]
      )
      # latent state predators, season 1
      pred_z[site,ds, 1] ~ dbern(
        pred_psi[site,ds,1]
      )
      for(t in 2:nseason){
        # logit linear predictor for predators, season >1
        logit(pred_psi[site,ds,t]) <- inprod(
          pred_beta[ds,],
          pred_psi_design_matrix[site,]
        ) + pred_theta[ds] * pred_z[site,ds, t-1]
        pred_z[site,ds,t] ~ dbern(
          pred_psi[site,ds,t]
        )
      }
    }
  }
  # Do the inxs stuff here as it is a lot of indexing
  for(ps in 1:nprey){
    for(ds in 1:npred){
      for(site in 1:nsite){
        for(tt in 1:nseason){
          for(i in 2:nparm_prey_psi){
            inxs_slopes[site,ps,ds,tt,i-1] <-  prey_psi_design_matrix[site,i] * pred_z[site,ds,tt] * inxs_beta[ps,ds,i-1]
          } 
        }
      }
    }
  }
  for(site in 1:nsite){
    for(ps in 1:nprey){
      logit(prey_psi[site, ps,1]) <- inprod(
        prey_beta[ps,],
        prey_psi_design_matrix[site, ]
      ) +
        inprod(
          inxs_b0[ps,],
          pred_z[site,,1]
        ) + sum(inxs_slopes[site,ps,1:npred,1,])
      # latent state prey, season 1
      prey_z[site,ps, 1] ~ dbern(
        prey_psi[site,ps,1]
      )
      for(t in 2:nseason){
        logit(prey_psi[site,ps,t]) <- inprod(
          prey_beta[ps,],
          prey_psi_design_matrix[site, ]
        ) +
          inprod(
            inxs_b0[ps,],
            pred_z[site,,1]
          ) + sum(inxs_slopes[site,ps,1:npred,1,]) +
          prey_theta[ps] * prey_z[site,ps,t-1]
        prey_z[site,ps,t] ~ dbern(
          prey_psi[site,ps,t]
        )
      }
    }
  }
  # Observational model predator
  for(ds in 1:npred){
    for(site in 1:nsite){
      for(tt in 1:nseason){
        logit(pred_rho[site,ds,tt]) <- inprod(
          pred_alpha[ds,],
          pred_rho_design_matrix[site,]
        )
        pred_y[site,ds, tt] ~ dbin(
          pred_rho[site,ds,tt] * pred_z[site,ds,tt],
          J[site,tt]
        )
      }
    }
  }
  # Observational model prey
  for(site in 1:nsite){
    for(ps in 1:nprey){
      for(tt in 1:nseason){
        logit(prey_rho[site,ps,tt]) <- inprod(
          prey_alpha[ps,],
          prey_rho_design_matrix[site,]
        )
        prey_y[site,ps, tt] ~ dbin(
          prey_rho[site,ps,tt] * prey_z[site,ps,tt],
          J[site,tt]
        )
      }
    }
  }
  # Priors
  for(ds in 1:npred){
    pred_theta[ds] ~ dlogis(0,1)
    for(i_pred_psi in 1:nparm_pred_psi){
      pred_beta[ds,i_pred_psi] ~ dlogis(0,1)
    }
    for(i_pred_rho in 1:nparm_pred_rho){
      pred_alpha[ds,i_pred_rho] ~ dlogis(0,1)
    }
    for(ps in 1:nprey){
      inxs_b0[ps,ds] ~ dlogis(0,1)
      for(iinxs_psi in 1:(nparm_prey_psi - 1)){
        inxs_beta[ps,ds,iinxs_psi] ~ dlogis(0,1)
      }
    }
  }
  for(ps in 1:nprey){
    prey_theta[ps] ~ dlogis(0,1)
    for(i_prey_psi in 1:nparm_prey_psi){
      prey_beta[ps, i_prey_psi] ~ dlogis(0,1)
    }
    for(i_prey_rho in 1:nparm_prey_rho){
      prey_alpha[ps, i_prey_rho] ~ dlogis(0,1)
    }
  }
}
# var prey_y[nsite,nprey, nseason], prey_rho[nsite,nprey,nseason], pred_y[nsite,npred, nseason], pred_rho[nsite,npred,nseason];
