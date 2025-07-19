## equations

# Force of infection equations
λ_fsw = c_fsw * β * (1 - ε_fsw * (1 - v)) * ((I_acute_clients + I_chronic_clients + I_chronic_ART_clients) / N_clients)

λ_fsw_prep = c_fsw * β * ο * (1 - ε_fsw * (1 - v)) * ((I_acute_clients + I_chronic_clients + I_chronic_ART_clients) / N_clients)

# Differential equations
dS_fsw_dt = γ_fsw - λ_fsw * S_fsw - PrEP_FSW * S_fsw - τ_fsw * S_fsw - d_fsw * S_fsw

dS_fsw_prep_dt = PrEP_FSW * S_fsw - λ_fsw_prep * S_fsw_prep - τ_fsw_prep * S_fsw_prep - d_fsw_prep * S_fsw_prep

dI_fsw_acute_dt = λ_fsw * S_fsw + λ_fsw_prep * S_fsw_prep - ζ_fsw_acute * I_fsw_acute - ζ_fsw_chronic_ART * I_fsw_acute - τ_fsw * I_acute_fsw - d_fsw * I_acute_fsw

dI_fsw_chronic_dt = ζ_fsw_acute * I_fsw_acute - ϕ_fsw_chronic * I_fsw_chronic - ν_fsw * I_fsw_chronic - τ_fsw * I_fsw_chronic - d_fsw * I_fsw_chronic

dI_fsw_chronic_ART_dt = (1 - ζ_fsw_acute) * I_fsw_acute - ϕ_fsw_chronic_ART * I_fsw_chronic_ART - ν_fsw_ART * I_fsw_chronic_ART - τ_fsw * I_fsw_chronic_ART - d_fsw * I_fsw_chronic_ART

dAIDS_dt = ν_fsw * I_fsw_chronic + ν_fsw_ART * I_fsw_chronic_ART - κ_fsw_AIDS * AIDS - τ_fsw * AIDS - d_fsw