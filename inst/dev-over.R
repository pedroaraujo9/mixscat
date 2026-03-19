library(tidyverse)
devtools::load_all()

ari = mclust::adjustedRandIndex

compute_asw = function(w, z_dist) {
  asw = cluster::silhouette(w, dist = z_dist)
  mean(asw[, 3])
}

get_w_ward_mode = function(M, model_data, stan_model, lambda) {

  w_ward_init = hclust(as.dist(model_data$data$z_dist), method = "ward.D") |> cutree(k = M)

  post_map = find_post_map(
    M = M,
    w = w_ward_init,
    stan_model = stan_model,
    init_list = NULL,
    lambda = lambda,
    model_data = model_data

  )

  list(
    w = w_ward_init,
    beta = cbind(post_map$par$beta, 0),
    pw = post_map$par$pw %>% as.numeric()
  )

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

M = 20
n_basis = 10
iters = 100
burn_in = 50
thin = 2
chains = 2
seed = 2

lambda = 1
intercept_penalty = 1
dirichlet_param = 0.001

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

devtools::load_all()


fits = readRDS("~/Documents/GitHub/clustering-mortality-two-step/models/mixscat-fit/mixscat-male-fit-wpp-G=5-lambda=1.rds")
w_init = fits$clusters[, 11]

fits = lapply(1:20, function(chain){

  print(paste("Chain", chain, "\n"))

  w_ward_init = get_w_ward_mode(
    M = chain, model_data = model_data, stan_model = stan_model, lambda = lambda
  )

  init_run = single_run(
    M = chain,
    w = NULL,
    model_data = model_data,
    lambda = 1,
    init_list = w_ward_init,
    iters = 100,
    burn_in = 0,
    thin = 1,
    verbose = TRUE,
    seed = NULL
  )

  w_init = init_run$sample_list$w[100, ]

  post_map = find_post_map(
    M = M,
    w = w_init,
    stan_model = stan_model,
    init_list = NULL,
    lambda = lambda,
    model_data = model_data
  )

  init_list = list(
    w_init = w_init,
    beta_init = cbind(post_map$par$beta, 0)
  )

  fit = single_run(
    M = M,
    w = NULL,
    model_data = model_data,
    lambda = 1,
    init_list = init_list,
    iters = 1000,
    burn_in = 0,
    thin = 1,
    verbose = TRUE,
    seed = NULL
  )

  fit

})


model_data$spline$dirichlet_param

size = lapply(1:length(fits), function(i){
  fits[[i]]$sample_list$w[iters,] %>% unique() %>% length()
}) %>% do.call(c, .)

logp_mode = lapply(1:length(fits), function(i){
  w = fits[[i]]$sample_list$w[1000,]

  post_map = find_post_map(
    M = M,
    w = w,
    stan_model = stan_model,
    init_list = NULL,
    lambda = lambda,
    model_data = model_data
  )

  post_map$par$log_like

})

llp = lapply(1:20, function(i){
  fits[[i]]$sample_list$logpost[500:1000, "logpost"] %>% mean()
}) %>% do.call(c, .)

plot(size, llp)

lapply(1:20, function(i){
  compute_asw(w = fits[[i]]$sample_list$w[1000, ], z_dist = model_data$data$z_dist)
}) %>%
  do.call(c, .) %>%
  plot()


fits[[9]]$sample_list$w[1000, ] %>% sort()


ari(compute_asw(w = fits[[i]]$sample_list$w[1000, ])


logp = lapply(1:length(fits), function(i){
  fits[[i]]$sample_list$logpost[iters, "logpost"]
}) %>% do.call(c, .)

lapply(1:length(fits), function(i){
  fits[[i]]$sample_list$logpost[500:1000, "logpost"] %>% mean()
}) %>% do.call(c, .) %>% plot()


fits[[1]]$sample_list$w[iters, ] %>% sort()

size %>% table() %>% prop.table()

fits[[1]]$sample_list$beta[, 1, 1] %>% plot()

w = fits[[1]]$sample_list$w[iters, ]

ari(w, fits[[1]]$init$w)

fits[[1]]$sample_list$logpost %>% colMeans()

lapply(1:1000, function(i){
  fits[[20]]$sample_list$w[i, ] %>% unique() %>% length()
}) %>% do.call(c, .) %>% plot(type="l")

lapply(100:1000, function(i){
  fits[[1]]$sample_list$logpost[i, "logpost"]
}) %>% do.call(c, .) %>% plot(type="l")

lapply(10:iters, function(i){
  fits[[1]]$sample_list$w[i, ] %>% compute_asw(model_data$data$z_dist)
}) %>% do.call(c, .) %>% plot()

fits[[1]]$sample_list$w[iters, ] %>% sort()


size %>% table() %>% plot()

w1 = fits[[3]]$sample_list$w[iters,]
w2 = fits[[4]]$sample_list$w[iters,]
w3 = fits[[9]]$sample_list$w[iters,]

ari(w1, w2)
ari(w1, w3)
ari(w2, w3)


fits[[3]]$sample_list$w[iters,] %>% unique() %>% length()

lapply(1:iters, function(i){
  fit$sample_list$w[i,] %>% unique() %>% length()
}) %>% do.call(c, .) %>% plot()


fit$sample_list$w[iters,] %>% unique() %>% length()
fit$sample_list$beta %>% compute_post_stat() %>% round(5)





