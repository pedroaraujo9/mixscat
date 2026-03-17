library(tidyverse)
devtools::load_all()

ari = mclust::adjustedRandIndex
compute_asw = function(w, z_dist) {
  asw = cluster::silhouette(w, dist = z_dist)
  mean(asw[, 3])
}

#### data ####
data = readRDS("~/Documents/GitHub/clustering-mortality-two-step/data/input-data-cluster.rds")

x_male = data %>%
  filter(sex == "Male") %>%
  select(-country, -year, -sex, -(Z2:Z6)) %>%
  as.matrix()

id = data %>% select(country, year) %>% distinct() %>% .$country
time = data %>% select(country, year) %>% distinct() %>% .$year

z_male = data %>%
  filter(sex == "Male") %>%
  as.data.frame() %>%
  .[, paste0("Z", 5)]

male_levels = z_male %>%
  unique() %>%
  sort() %>%
  as.character()

z = factor(z_male, levels = male_levels)
n_id = length(unique(id))
n_time = length(unique(time))
id_unique = id |> unique() |> sort()

z_seq = z |>
  matrix(nrow = n_id, ncol = n_time, byrow = T) |>
  TraMineR::seqdef()

rownames(z_seq) = id_unique
z_dist = TraMineR::seqdist(z_seq, method = "HAM")


data.frame(
  id = id,
  time = time,
  z = z
) %>%
  ggplot(aes(x=time, y=id, fill=z)) +
  geom_tile(color="black") +
  viridis::scale_fill_viridis(discrete = T)

#### ward ####
get_w_ward = function(M, z_dist) {
  w_ward_init = hclust(as.dist(z_dist), method = "ward.D") |> cutree(k = M)
  w_ward_init
}

ward_asw = lapply(2:20, function(m){
  w_ward_init = get_w_ward(m, z_dist)
  compute_asw(w_ward_init, z_dist)
}) %>% do.call(c, .)

#### MED seq ####
get_w_med = function(M, z_seq) {

  w_med_init = MEDseq::MEDseq_fit(z_seq, G = M, modtype = "UU") |>
    MEDseq::get_MEDseq_results(what = "MAP")

  return(w_med_init)
}

med_asw = lapply(2:20, function(m){
  w_med_init = get_w_med(m, z_seq)
  compute_asw(w_med_init, z_dist)
}) %>% do.call(c, .)

data.frame(
  med_asw = med_asw,
  ward_asw = ward_asw,
  M = 2:20
) %>%
  gather(method, asw, -M) %>%
  ggplot(aes(x=M, y=asw, color=method)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = 2:20)

data %>%
  filter(sex == "Male", country %in% c("Argentina","Bulgaria","Malaysia","Serbia","Uruguay")) %>%
  ggplot(aes(x=year, y=country, fill=factor(Z5))) +
  geom_tile(color="black")


#### data model ####
M = 11
n_basis = 10
iters = 100
burn_in = 50
thin = 2
chains = 2
seed = 2

lambda = 5
intercept_penalty = 1
dirichlet_param = 1

extraDistr::rdirichlet(1000, rep(0.05, M)) %>% min()

model_data = create_model_data(
  z = z,
  id = id,
  time = time,
  n_basis = n_basis,
  intercept_penalty = intercept_penalty,
  dirichlet_param = dirichlet_param
)


model_path = system.file(
  "stan", "stan_full_logpost.rds", package = "mixscat"
)

stan_model = readRDS(model_path)

#### inits ####
w_ward_init = get_w_ward(M, z_dist)

w_init_opt = find_post_map(
  M = M,
  w = w_ward_init,
  stan_model = stan_model,
  init_list = NULL,
  lambda = lambda,
  model_data = model_data
)

beta_init_ward = w_init_opt$par$beta
pw_init_ward = w_init_opt$par$pw

iters = 200
ward_init_runs = lapply(1:5, function(i){

  cat("Init run:", i, "\r")

  run = single_run(
    M = M,
    w = NULL,
    model_data = model_data,
    lambda = lambda,
    init_list = list(
      beta = cbind(beta_init_ward, 0),
      pw = pw_init_ward,
      w = w_ward_init
    ),
    iters = iters,
    burn_in = 0,
    thin = 1,
    verbose = FALSE,
    seed = NULL
  )

  w = run$sample_list$w[iters, ]

  opt_run = find_post_map(
    M = M,
    w = w,
    stan_model = stan_model,
    init_list = NULL,
    lambda = lambda,
    model_data = model_data
  )

  logpost = opt_run$par$log_like

  list(
    logpost = logpost,
    beta = opt_run$par$beta,
    pw = opt_run$par$pw,
    w = w
  )

})

ward_init_runs %>% purrr::map_dbl(~.x$logpost)

mclust::adjustedRandIndex(
  ward_init_runs[[1]]$w,
  ward_init_runs[[2]]$w
)


idx = ward_init_runs %>% purrr::map_dbl(~.x$logpost) %>% which.max()
w_ward_random_init = ward_init_runs[[idx]]$w
beta_ward_random_init = ward_init_runs[[idx]]$beta
pw_ward_random_init = ward_init_runs[[idx]]$pw

sd_beta = compute_beta_sd_matrix(model_data, lambda, M = M)

compute_logpost(
  beta = cbind(beta_ward_random_init, 0),
  pw = as.numeric(pw_ward_random_init),
  w = w_ward_random_init,
  model_data = model_data,
  sd_beta = sd_beta
)

compute_logpost(
  beta = cbind(beta_init_ward, 0),
  pw = as.numeric(pw_init_ward),
  w = w_ward_init,
  model_data = model_data,
  sd_beta = sd_beta
)

ari(w_ward_random_init, w_ward_init)
table(w_ward_random_init, w_ward_init)


init_list = find_init_w(
  M = 20,
  model_data = model_data,
  n_init = 10,
  init_mcmc_iters = 100,
  lambda = 5,
  seed = 1,
  verbose = TRUE
)

init_list$beta_init
init_list$ari_init

M = 2
chains = 3
iters = 200
burn_in = 100
thin = 2
lambda = 1
dirichlet_param = 1
intercept_penalty = 1
seed = NULL
n_init = 10
init_mcmc_iters = 50
verbose = TRUE
n_basis = 10

model_data = create_model_data(
  z = z,
  id = id,
  time = time,
  n_basis = n_basis,
  intercept_penalty = intercept_penalty,
  dirichlet_param = dirichlet_param
)

run = pipeline(
  model_data = model_data,
  M = M,
  chains = chains,
  iters = iters,
  burn_in = burn_in,
  thin = thin,
  lambda = 5,
  dirichlet_param = dirichlet_param,
  seed = seed,
  n_init = n_init,
  init_mcmc_iters = init_mcmc_iters,
  verbose = TRUE
)

run$metrics
run$conv_check$logpost

M = 2
chains = 2
iters = 20
burn_in = 10
thin = 2
lambda = 1
dirichlet_param = 1
intercept_penalty = 1
seed = 1
n_init = 2
init_mcmc_iters = 50
verbose = TRUE
n_basis = 10
n_cores = 2


de
fit = fit_mixscat(
  M = 1,
  z = z,
  id = id,
  time = time,
  n_basis = n_basis,
  iters = iters,
  thin = thin,
  burn_in = burn_in,
  chains = chains,
  lambda = lambda,
  intercept_penalty = intercept_penalty,
  dirichlet_param = dirichlet_param,
  n_init = n_init,
  init_mcmc_iters = init_mcmc_iters,
  n_cores = n_cores,
  seed = seed,
  verbose = TRUE
)


runs = lapply(1:2, function(i){
  single_run(
    M = 11,
    w = NULL,
    model_data = model_data,
    lambda = lambda,
    init_list = init_list,
    iters = 500,
    burn_in = 100,
    thin = 2,
    verbose = TRUE,
    seed = NULL
  )
})

w1 = runs[[1]]$sample_list$w %>% comp_class()
w2 = runs[[2]]$sample_list$w %>% comp_class()
ari(w1, w2)

ari(init_list$w_init, w1)
ari(init_list$w_init, w2)

runs[[1]]$sample_list$logpost[, "loglike"] %>% plot()
runs[[2]]$sample_list$logpost[, "loglike"] %>% plot()

data.frame(
  ll1 = runs[[1]]$sample_list$logpost[, "loglike"],
  ll2 = runs[[2]]$sample_list$logpost[, "loglike"]
) %>%
  mutate(iter = 1:nrow(.)) %>%
  gather(run, loglike, -iter) %>%
  ggplot(aes(x=iter, y=loglike, color=run)) +
  geom_line()

posterior::rhat(rbind(
  runs[[1]]$sample_list$logpost[, "loglike"],
  runs[[2]]$sample_list$logpost[, "loglike"]
))

posterior::rhat(rbind(
  runs[[1]]$sample_list$logpost[, "logpost"],
  runs[[2]]$sample_list$logpost[, "logpost"]
))

runs[[1]]$sample_list$logpost[, "logpost"] %>% mean()
runs[[2]]$sample_list$logpost[, "logpost"] %>% mean()

ari(w1, w2)
compute_asw(w1, z_dist)
compute_asw(w_ward_init, z_dist)

table(w1, w2)

ms = 1:20
lambda = 1

devtools::load_all()

model_path = system.file(
  "stan", "stan_full_logpost.rds", package = "mixscat"
)

stan_model = readRDS(model_path)

w_init_mode = lapply(ms, function(m){

  cat("Finding mode for m =", m, "\r")

  opt_run = find_post_map(
    M = 20,
    w = get_w_ward(m, z_dist),
    stan_model = stan_model,
    init_list = NULL,
    lambda = lambda,
    model_data = model_data
  )

  opt_mode = opt_run$par$log_like

  list(
    opt_mode = opt_mode # , mcmc_mode = mcmc_mode
  )

})

logpost = w_init_mode %>% purrr::map_dbl(~{.x$opt_mode})
plot(logpost)
which.max(logpost)

run_init = single_run(
  M = M,
  w = w_ward_init,
  model_data = model_data,
  lambda = lambda,
  init_list = NULL,
  iters = 200,
  burn_in = 100,
  thin = 10,
  verbose = TRUE,
  seed = NULL
)

run_init$sample_list$beta

beta_init = run_init$sample_list$beta %>% compute_post_stat()
pw_init = run_init$sample_list$pw %>% compute_post_stat()

init_list = list(
  beta = beta_init,
  pw = pw_init,
  w = w_init_matrix[, 11]
)




w_init_matrix = lapply(1:20, function(m){
  hclust(as.dist(model_data$data$z_dist), method = "ward.D") |> cutree(m)
}) %>% do.call(cbind, .)

apply(w_init_matrix[,-1], 2, function(w) compute_asw(w, model_data)) %>% plot()

lambda_opt = lapply(1:20, function(m){

  cat("Calibrating lambda for m =", m, "\r")

  lambda_opt = calibrate_lambda(
    w = hclust(as.dist(model_data$data$z_dist)) |> cutree(m),
    M = m,
    model_data = model_data,
    lambda_grid = seq(from = 0.01, to = 5, length.out = 30)
  )

  lambda_opt$best_lambda

}) %>% do.call(c, .)

lambda_opt

lambda = lambda_opt$best_lambda

model_path = system.file(
  "stan", "stan_full_logpost.rds", package = "mixscat"
)

stan_model = readRDS(model_path)

lambda = 1

ms = 1:20
devtools::load_all()
w_init_mode = lapply(ms, function(m){

  cat("Finding mode for m =", m, "\r")

  opt_run = find_post_map(
    M = M,
    w = w_init[, m],
    stan_model = stan_model,
    init_list = NULL,
    lambda = lambda,
    model_data = model_data
  )

  opt_mode = opt_run$par$log_like

  # mcmc_run = single_run(
  #   M = M,
  #   w = w_init[, m],
  #   model_data = model_data,
  #   lambda = lambda,
  #   init_list = NULL,
  #   iters = 100,
  #   burn_in = 50,
  #   thin = 2,
  #   verbose = TRUE,
  #   seed = 1
  # )
  #
  # mcmc_mode = mcmc_run$sample_list$logpost[, 2] %>% mean()

  list(
    opt_mode = opt_mode # , mcmc_mode = mcmc_mode
  )

})

#mcmc_mode = purrr::map_dbl(w_init_mode, ~.x$mcmc_mode)
opt_mode = purrr::map_dbl(w_init_mode, ~.x$opt_mode)
#plot(10:20, mcmc_mode)
plot(ms[5:20], opt_mode[5:20])

which.max(mcmc_mode)
which.max(opt_mode)

log_post = purrr::map_dbl(w_init_mode, ~.x$par$log_like)
plot(log_post)


which.max(log_post)

fit$par$log_like


run_init = single_run(
  M = M,
  w = w_init_matrix[, 11],
  model_data = model_data,
  lambda = lambda,
  init_list = NULL,
  iters = 200,
  burn_in = 100,
  thin = 10,
  verbose = TRUE,
  seed = NULL
)

beta_init = run_init$sample_list$beta %>% compute_post_stat()
pw_init = run_init$sample_list$pw %>% compute_post_stat()

init_list = list(
  beta = beta_init,
  pw = pw_init,
  w = w_init_matrix[, 11]
)

# init_list = list(
#   w = rep(1, model_data$dims$n_id)
# )

iters = 300
cut_point = 1
temp_vec = c(
  seq(1, 1, length.out= cut_point),
  rep(1, iters - cut_point)
)
plot(temp_vec)

devtools::load_all()
runs = lapply(1:3, function(i) {
  if(i == 1) {
    single_run(
      M = M,
      w = NULL,
      model_data = model_data,
      lambda = lambda,
      init_list = NULL,
      iters = iters,
      burn_in = 0,
      thin = 1,
      verbose = TRUE,
      seed = NULL,
      temperature_vec = temp_vec
    )
  }else{
    single_run(
      M = M,
      w = NULL,
      model_data = model_data,
      lambda = lambda,
      init_list = init_list,
      iters = iters,
      burn_in = 0,
      thin = 1,
      verbose = TRUE,
      seed = NULL,
      temperature_vec = temp_vec
    )
  }
})

w1 = runs[[1]]$sample_list$w[iters, ]
w2 = runs[[2]]$sample_list$w[iters, ]
w3 = runs[[3]]$sample_list$w[iters, ]

runs[[1]]$sample_list$logpost[iters, ]
runs[[2]]$sample_list$logpost[iters, ]
runs[[3]]$sample_list$logpost[iters, ]

ari = mclust::adjustedRandIndex

ari(w1, w2)
table(w1, w2)

ari(w1, w_init_matrix[, 11])
ari(w2, w_init_matrix[, 11])
ari(w3, w_init_matrix[, 11])

w_med = model_data$data$z |>
  matrix(nrow = model_data$dims$n_id, ncol = model_data$dims$n_time, byrow = T) |>
  TraMineR::seqdef() %>%
  MEDseq::MEDseq_fit(G = M, modtype = "CC") |>
  MEDseq::get_MEDseq_results(what = "MAP")

compute_asw(w1, model_data)
compute_asw(w2, model_data)
compute_asw(w_init_matrix[,11], model_data)
compute_asw(w_med, model_data)

names(w_med) = names(w1)
w_med %>% sort() %>% unique()

init1 = find_init_w(
  M = M,
  model_data = model_data,
  seed = 1,
  init_control = init_control
)

init2 = find_init_w(
  M = M,
  model_data = model_data,
  seed = 2,
  init_control = init_control
)

init1$logpost
init2$logpost

w1 = init1$w
w2 = init2$w

ari(w1, w2)
w1 %>% sort()
w2 %>% sort()

table(w1, w2)

devtools::load_all()
init_runs = lapply(1, function(seed) {
  find_init_w(
    M = M,
    model_data = model_data,
    seed = seed,
    init_control = init_control
  )
})


init_runs[[1]]$logpost
init_runs[[1]]$w %>% unique() %>% length()
init_runs[[1]]$w %>% sort()

init_runs[[1]]$sample_list$w %>% comp_class() %>% sort()


init_runs[[1]]$logpost
init_runs[[2]]$logpost
init_runs[[3]]$logpost

init_runs[[1]]$w %>% unique() %>% length()
init_runs[[2]]$w %>% unique() %>% length()
init_runs[[3]]$w %>% unique() %>% length()

mclust::adjustedRandIndex(init_runs[[1]]$w, init_runs[[2]]$w)
mclust::adjustedRandIndex(init_runs[[1]]$w, init_runs[[3]]$w)
mclust::adjustedRandIndex(init_runs[[2]]$w, init_runs[[3]]$w)

init_runs[[1]]$w %>% sort()

