library(tidyverse)
devtools::load_all()

ari = mclust::adjustedRandIndex
compute_asw = function(w, z_dist) {
  asw = cluster::silhouette(w, dist = z_dist)
  mean(asw[, 3])
}

#### data ####
data = readRDS("~/Documents/GitHub/clustering-mortality-two-step/data/input-data-cluster.rds")
sex_str = "Male"
G = 5

x = data %>%
  filter(sex == sex_str) %>%
  select(-country, -year, -sex, -(Z2:Z6)) %>%
  as.matrix()

id = data %>% select(country, year) %>% distinct() %>% .$country
time = data %>% select(country, year) %>% distinct() %>% .$year

z = data %>%
  filter(sex == sex_str) %>%
  as.data.frame() %>%
  .[, paste0("Z", G)]

levels = z %>%
  unique() %>%
  sort() %>%
  as.character()

z = factor(z, levels = levels)
n_id = length(unique(id))
n_time = length(unique(time))
id_unique = id |> unique() |> sort()

data.frame(
  id = id,
  time = time,
  z = z
) %>%
  ggplot(aes(x=time, y=id, fill=z)) +
  geom_tile(color="black") +
  viridis::scale_fill_viridis(discrete = T)


data %>%
  filter(sex == "Male", country %in% c("Argentina","Bulgaria","Malaysia","Serbia","Uruguay")) %>%
  ggplot(aes(x=year, y=country, fill=factor(Z5))) +
  geom_tile(color="black")


#### data model ####
M = 15
n_basis = 10
iters = 100
burn_in = 50
thin = 2
chains = 2
seed = 1


n_cores = 1
verbose = TRUE

lambda = 1
intercept_penalty = 1
dirichlet_param = 0.01
submodels = 1:M

model_data = create_model_data(
  z = z,
  id = id,
  time = time,
  n_basis = n_basis,
  intercept_penalty = intercept_penalty,
  dirichlet_param = dirichlet_param
)

iters = 20

devtools::load_all()
n_clust = find_number_clust(
  M_max = 15,
  n_init = 50,
  iters = iters,
  lambda_grid = c(1),
  dirichlet_param = 0.01,
  intercept_penalty = 1,
  model_data,
  seed = 1,
  verbose = TRUE
)


n_clust_distri = lapply(names(n_clust), function(lambda){

  p = n_clust[[lambda]]$n_clust_chain[, iters] %>% table() %>% prop.table()

  data.frame(
    n_clust = as.numeric(names(p)),
    prop = as.numeric(p),
    lambda = lambda
  )

}) %>% do.call(rbind, .)

n_clust_distri %>%
  ggplot(aes(x=n_clust, y=prop, group=lambda, color=factor(lambda, levels = c(1, 2, 3, 5)))) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = 1:20)


total_distri = lapply(names(n_clust), function(lambda){

  iters = n_clust[[lambda]]$n_clust_chain %>% ncol()
  n_clust[[lambda]]$n_clust_chain[, iters]

}) %>% do.call(c, .) %>% table() %>% prop.table()

barplot(total_distri)

w_sample = lapply(names(n_clust), function(lambda){

  n_clust[[lambda]]$w

}) %>% do.call(rbind, .)

w_post_prob = lapply(names(n_clust), function(lambda){

  n_clust[[lambda]]$post_prob

})

psm = mcclust::comp.psm(cls = w_sample)
cl = mcclust::maxpear(psm, cls.draw = w_sample, max.k = 15, method = "avg")
colnames(psm) = rownames(psm) = model_data$data$id_unique
w_pear = cl$cl
names(w_pear) = model_data$data$id_unique
w_pear %>% unique() %>% length()

w_sample


ls = label.switching::label.switching(
  method = "STEPHENS",
  z = w_sample,
  p = w_post_prob[[1]],
  K = 15
)

new_w_sample = w_sample
for(i in 1:nrow(w_sample)){
  new_w_sample[i, ] = ls$permutations$STEPHENS[i, w_sample[i, ]]
}

lapply(1:nrow(new_w_sample), function(i){
  new_w_sample[i, ] %>% unique() %>% length()
}) %>% do.call(c, .) %>% table() %>% barplot()

new_w_sample[, "Germany"] %>% table() %>% prop.table()


cl = mcclust::maxpear(psm, cls.draw = w_sample, max.k = 15, method = "draws")
colnames(psm) = rownames(psm) = model_data$data$id_unique
w_pear = cl$cl
names(w_pear) = model_data$data$id_unique
w_pear %>% sort() %>% unique() %>% length()

w_pear

w_pear %>% sort()
w_pear %>% unique() %>% length()

ls = label.switching::label.switching(
  method = "ECR",
  z = w_sample,
  zpivot = w_pear,
  K = 15
)

permu = ls$permutations$ECR

lapply(1:120, function(iter){
  w_sample[iter, ] %>% unique() %>% length()
}) %>% do.call(c, .) %>% table() %>% barplot()
w_sample %>% dim()



w_sample_new = w_sample
for(i in 1:nrow(w_sample)){
  w_sample_new[i, ] = permu[i, w_sample[i, ]]
}

colnames(w_sample_new) = model_data$data$id_unique

w_sample_new[, "Argentina"] %>% table() %>% prop.table() %>% plot()

post_prob = lapply(1:ncol(w_sample_new), function(i){
  w_sample_new[, i] %>% table() %>% prop.table() %>% as.numeric()
}) %>% do.call(rbind, .)

rownames(post_prob) = model_data$data$id_unique
View(post_prob)

post_prob[, 14] %>% sort()

w = extraDistr::rcat(n = nrow(post_prob), prob = post_prob)

names(w) = model_data$data$id_unique

w %>% sort()


w_sample_new[, 1] %>% table() %>% prop.table()

mean(w_pear  == ls$clusters[1,])

cl = ls$clusters %>% as.numeric()
names(cl) = model_data$data$id_unique
cl %>% sort()

colnames(wf) = model_data$data$id_unique
mean(wf[, "Uruguay"] == wf[, "Bulgaria"])

cl %>% sort()

iters = 1000
burn_in = 0
thin = 1
n_runs = 1
lambda = 10
sqrt(1/lambda)

devtools::load_all()
runs = lapply(1:n_runs, function(i){

  single_run(
    M = M,
    w = NULL,
    model_data = model_data,
    lambda = lambda,
    intercept_penalty = intercept_penalty,
    dirichlet_param = dirichlet_param,
    init_list = list(w = get_w_ward(M = M, z_dist = model_data$data$z_dist)),
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    verbose = verbose,
    add_logpost = TRUE,
    seed = NULL
  )

})

runs[[1]]$sample_list$beta[-c(1:100), 1, 1] %>% plot()
beta = runs[[1]]$sample_list$beta[floor(iters/2):iters,,] %>% compute_post_stat() %>% round(5)
beta[1, ]
mclust::softmax(beta[1, ])
model_data$data$z %>% table() %>% prop.table()
beta



lapply(1:iters, function(i){
  runs[[1]]$sample_list$w[i, ] %>% unique() %>% length()
}) %>% do.call(c, .) %>% plot(type = "l")

lapply(1:iters, function(i){
  runs[[2]]$sample_list$w[i, ] %>% unique() %>% length()
}) %>% do.call(c, .) %>% plot(type = "l")

lapply(1:iters, function(i){
  runs[[3]]$sample_list$w[i, ] %>% unique() %>% length()
}) %>% do.call(c, .) %>% plot(type = "l")


w1 = runs[[1]]$sample_list$w[floor(iters/2):iters,] %>% comp_class()
w2 = runs[[2]]$sample_list$w[floor(iters/2):iters,] %>% comp_class()
w3 = runs[[3]]$sample_list$w[floor(iters/2):iters,] %>% comp_class()

ari(w1, w2)
ari(w1, w3)
ari(w2, w3)

table(w1, w2)
names(w3[w3 == 5])[names(w3[w3 == 5]) %in% names(w2[w2 == 4])]

runs[[3]]$sample_list$w[,"Colombia"] %>% table() %>% prop.table()
runs[[2]]$sample_list$w[,"Colombia"] %>% table() %>% prop.table()

lapply(1:iters, function(i){
  runs[[2]]$sample_list$w[i, ] %>% unique() %>% length()
}) %>% do.call(c, .) %>% plot(type = "l")


loglike_bayes = compute_logpost(
  beta = runs[[1]]$sample_list$beta[floor(iters/2):iters,,] %>% compute_post_stat(),
  w = runs[[1]]$sample_list$w[floor(iters/2):iters,] %>% comp_class(),
  pw = runs[[1]]$sample_list$pw[floor(iters/2):iters,] %>% compute_post_stat(),
  model_data = model_data,
  dirichlet_param = dirichlet_param,
  beta_sd = runs[[1]]$beta_sd
)["loglike"]

loglike_mean = mean(runs[[1]]$sample_list$logpost[floor(iters/2):iters, "loglike"])

p = 2 * (loglike_bayes - loglike_mean)
p


l = runs[[1]]$sample_list$logpost[200:300, "loglike"]
-2*mean(l) + log(71) * p
mean(l)


table(w1, w2)

beta = run$sample_list$beta[floor(iters/2):iters,,] %>% compute_post_stat() %>% round(5)
w = run$sample_list$w[floor(iters/2):iters,] %>% comp_class()
run$sample_list$w_post_prob %>% compute_post_stat() %>% round(5)



beta[1, ] %>% mclust::softmax()
model_data$data$z %>% table() %>% prop.table()


logit = function(x) log(x/(1-x))
run$sample_list$w_post_prob %>% compute_post_stat()

B = cbind(1, model_data$spline$B)

#### runs over lambda ####


model_path = system.file(
  "stan", "stan_full_logpost.rds", package = "mixscat"
)

stan_model = readRDS(model_path)

iters = 100
burn_in = 0
thin = 1
n_runs = 1
lambda = 10
M = 15


runs = lapply(c(1, 2, 3, 4, 5), function(lambda){

  w = get_w_ward(M = M, z_dist = model_data$data$z_dist)

  post_map = find_post_map(
    M = M,
    w = get_w_ward(M = M, z_dist = model_data$data$z_dist),
    stan_model = stan_model,
    init_list = NULL,
    lambda = lambda,
    intercept_penalty = intercept_penalty,
    dirichlet_param = dirichlet_param,
    model_data = model_data
  )

  init_list = list(
    beta = cbind(post_map$par$beta, 0),
    w = w
  )


  single_run(
    M = M,
    w = NULL,
    model_data = model_data,
    lambda = lambda,
    intercept_penalty = intercept_penalty,
    dirichlet_param = dirichlet_param,
    init_list = init_list,
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    verbose = verbose,
    add_logpost = TRUE,
    seed = NULL
  )

})

lapply(1:length(runs), function(i){
  runs[[i]]$sample_list$w[iters,] %>% unique() %>% length()
}) %>% do.call(c, .)


lapply(1:iters, function(i){
  runs[[6]]$sample_list$w[i, ] %>% unique() %>% length()
}) %>% do.call(c, .)

runs = lapply(1:3, function(i){

  single_run(
    M = 6,
    w = NULL,
    model_data = model_data,
    lambda = 1,
    intercept_penalty = intercept_penalty,
    dirichlet_param = 1,
    init_list = list(w = get_w_ward(6, z_dist = model_data$data$z_dist)),
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    verbose = verbose,
    add_logpost = TRUE,
    seed = NULL
  )

})


w1 = runs[[1]]$sample_list$w %>% comp_class()
w2 = runs[[2]]$sample_list$w %>% comp_class()
w3 = runs[[3]]$sample_list$w %>% comp_class()
runs[[1]]$sample_list$logpost[-c(1:100), "loglike"] %>% plot()
runs[[2]]$sample_list$logpost[-c(1:100), "loglike"] %>% plot()

ari(w1, w2)
ari(w1, w3)
ari(w2, w3)

n_runs = 1
iters = 1000
burn_in = 0
thin = 1
lambda_grid = seq(100, 1, length.out = iters)
lambda_grid
#lambda_grid[floor(iters/2):iters] = 1
t = 1:iters
lambda_max = 100
lambda_min = 1
a = 0.01
t0 = iters / 3

rho = 1 / (1 + exp(-a * (t - t0)))

lambda = lambda_max^(1 - rho) * lambda_min^rho

lambda %>% tail()


lambda_grid = 1 + exp(3 - 0.01 * (1:iters))
lambda_grid %>% tail()
plot(lambda_grid)

n_runs = 20
iters = 20
lambda = 15
M = 15

runs = lapply(1:n_runs, function(i){

  single_run(
    M = M,
    w = NULL,
    model_data = model_data,
    lambda = lambda,
    intercept_penalty = intercept_penalty,
    dirichlet_param = dirichlet_param,
    init_list = NULL,
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    verbose = verbose,
    add_logpost = FALSE,
    seed = NULL
  )

})

lapply()

n_clust_iters = lapply(1:n_runs, function(i){

  lapply(1:iters, function(j){
    runs[[i]]$sample_list$w[j, ] %>% unique() %>% length()
  }) %>% do.call(c, .)

}) %>%
  do.call(rbind, .) %>%
  as.data.frame() %>%
  mutate(run = 1:nrow(.)) %>%
  gather(iter, n_clust, -run) %>%
  mutate(iter = as.numeric(gsub("V", "", iter))) %>%
  as_tibble()

n_clust_iters %>%
  filter(iter == max(iter)) %>%
  group_by(n_clust) %>%
  summarise(n = n()) %>%
  right_join(
    n_clust_iters %>%
      filter(iter == max(iter)),
    by = "n_clust"
  ) %>%
  select(run, n) %>%
  right_join(
    n_clust_iters,
    by = "run"
  ) %>%
  ggplot(aes(x=iter, y=n_clust, group=run, color=factor(n))) +
  geom_step(alpha = 1)


lapply(1:iters, function(i){
  runs[[4]]$sample_list$w[i, ] %>% unique() %>% length()
}) %>% do.call(c, .)

n_clust = lapply(1:iters, function(i){
  runs[[1]]$sample_list$w[i, ] %>% unique() %>% length()
}) %>% do.call(c, .)

n_clust

plot(n_clust)

plot(x = lambda_grid, y = n_clust)


