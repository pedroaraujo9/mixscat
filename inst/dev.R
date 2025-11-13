devtools::load_all()
devtools::document()

data = sim_data$data
time = data$time
id = data$id
x = NULL
w = NULL
z = sim_data$data$true_z
G = 3
M = 3
n_basis = 5
intercept = FALSE

model_data = create_model_data(
  time = time,
  id = id,
  x = x,
  z = z,
  w = w,
  G = G,
  M = M,
  n_basis = n_basis,
  intercept = intercept
)

priors = list(
  epsilon_w = 1,
  beta_sd = sqrt(10),
  mu_sd = sqrt(10),
  sigma_a = 1,
  sigma_b = 1
)

config = list(
  bounds = c(0.01, 10),
  n_points = 10,
  n_start_iters = NULL,
  epsilon_w = 1,
  beta_sd = sqrt(10),
  mu_sd = sqrt(10),
  sigma_a = 1,
  sigma_b = 1
)

#### single core ####
verbose = FALSE
seed = 2
lambda = 1
n_start = 5
iters = 20
n_cores = 1

inits = find_init(
  n_start = n_start,
  iters = iters,
  n_cores = n_cores,
  model_data = model_data,
  lambda = lambda,
  init_list = init_list,
  priors = priors,
  seed = seed
) |>
  expect_no_error() |>
  expect_no_message()


w = inits$init_list$w[1,]

z_ham_dist = compute_hamming(z = z, model_data = model_data)
w = hclust(z_ham_dist, method = "ward.D") |> cutree(k = model_data$M)

mclust::adjustedRandIndex(w, sim_data$true_w)
mclust::adjustedRandIndex(w2, sim_data$true_w)


lambda_grid = seq(config$bounds[1], config$bounds[2], length.out = config$n_points)
lambda_grid


obj_fun = function(lambda) {
  fit = find_map(
    z = z,
    w = w,
    lambda = lambda,
    stan_model = stan_model,
    model_data = model_data,
    fixed_sd = config$beta_sd
  )

  return(list(Score = fit$opt_fit$par$penal_ll))

}

penal = purrr::map_dbl(lambda_grid, ~{
  obj_fun(.x)$Score
})

best_lambda = lambda_grid[which.max(penal)]
best_lambda


G = 3
M = c(2, 3, 4)
z = sim_data$data$true_z %>% factor(labels = c("A", "B", "C"))
w = NULL
x = NULL
id = sim_data$data$id
time = sim_data$data$time
iters = 500
burn_in = iters/2
thin = 2
lambda = NULL
n_basis = 10
init_list = NULL
n_cores = 1
config = list(
  bounds = c(0.01, 5),
  lambda_start = 1,
  n_points = 20,
  n_start = 20,
  n_start_iters = 10,
  n_start_cores = 1,
  epsilon_w = 1,
  beta_sd = sqrt(10),
  mu_sd = sqrt(10),
  sigma_a = 1,
  sigma_b = 1,
  single_group = FALSE
)
verbose = TRUE
seed = NULL

fit = fit_mbfseq(
  G = G,
  M = M,
  z = z,
  w = w,
  x = x,
  id = id,
  time = time,
  iters = iters,
  burn_in = burn_in,
  thin = thin,
  lambda = lambda,
  n_basis = n_basis,
  init_list = init_list,
  n_cores = n_cores,
  config = config,
  verbose = verbose,
  seed = seed
)

fit$lambda_opt$`G=3, M=2`$lambda_plot
fit$lambda_opt$`G=3, M=3`$lambda_plot
fit$metrics$DIC3 %>% plot()

fit %>% plot_seq_cluster(M = 3, G = 3)

model_data = create_model_data(
  time = time,
  id = id,
  x = x,
  z = z,
  w = w,
  G = G,
  M = 3,
  n_basis = n_basis,
  intercept = FALSE
)



fit$metrics$AICM %>% plot()
fit$metrics$l_mean %>% plot()
2 * fit$metrics$l_var
10 * 2:4 * 2

D = -2 * fit$models$`G=3, M=3`$logpost[, 1]
D1 = mean(D)

alpha_post_mean = fit$models$`G=3, M=3`$sample_list$alpha %>% comp_post_stat(mean)
lambda = fit$lambda_opt$`G=3, M=3`$best_lambda
z = fit$models$`G=3, M=3`$z_class
w = fit$models$`G=3, M=3`$w_class
fixed_sd = config$beta_sd

model_data = create_model_data(
  time = time,
  id = id,
  x = x,
  z = z,
  w = w,
  G = G,
  M = 4,
  n_basis = n_basis,
  intercept = FALSE
)

D2 = -2 * eval_seq_logpost(
  z = z,
  alpha = alpha_post_mean,
  w = w, lambda = lambda, fixed_sd = fixed_sd, model_data = model_data
)

D1 - D2

model.frame.default()




post_summ = fit %>% posterior_summary(M = 3)
post_summ$`G=3, M=3` %>% plot_probability()

fit$model_data$z_levels


devtools::load_all()
fit %>%
  plot_seq_cluster(M = 2, G = 3,
                   cluster_label = c("Cluster\n  aa", "Cluster\nbb"),
                   label_hjust = -0.3)
fit %>% plot_seq_cluster(M = 3, G = 3)


#### real data ####
library(tidyverse)
library(mbfseq)
data = readRDS("~/Documents/GitHub/clustering-mortality-two-step/data/input-data-cluster.rds")

x_male = data %>%
  filter(sex == "Male") %>%
  select(-country, -year, -sex, -Z) %>%
  as.matrix()

x_female = data %>%
  filter(sex == "Female") %>%
  select(-country, -year, -Adult, -sex, -Z) %>%
  as.matrix()

id = data %>% select(country, year) %>% distinct() %>% .$country
time = data %>% select(country, year) %>% distinct() %>% .$year

z_male = data %>%
  filter(sex == "Male") %>%
  .$Z %>%
  factor(labels = c("M High", "M High Adult", "M Mid"))

z_female = data %>%
  filter(sex == "Female") %>%
  .$Z %>%
  factor(levels = c("F High", "F Mid", "F Low"))

z_female
#### model fit ####
config = list(
  single_group = FALSE,
  bounds = c(0.01,  5),
  lambda_start = 1,
  n_points = 20,
  n_start = 20,
  n_start_iters = 30,
  n_start_cores = 1,
  epsilon_w = 1,
  beta_sd = 1,
  mu_sd = 1,
  sigma_a = 1, sigma_b = 1
)

iters = 500
burn_in = 100
thin = 2
n_cores = 1
n_basis = 10
seed = 1
M = 2:6
lambda = NULL

male_fit = fit_mbfseq(
  M = M,
  z = z_male,
  id = id,
  w = NULL,
  x = NULL,
  time = time,
  iters = iters,
  burn_in = burn_in,
  thin = thin,
  n_cores = n_cores,
  n_basis = n_basis,
  seed = seed,
  verbose = TRUE,
  config = config,
  lambda = lambda
)

male_fit$metrics$DIC2 %>% plot()
male_fit$metrics$DIC3 %>% plot()



library(tidyverse)
devtools::load_all()
#### data ####
data = readRDS("../clustering-mortality-two-step/data/input-data-cluster.rds")

x_male = data %>%
  filter(sex == "Male") %>%
  select(-country, -year, -sex, -Z) %>%
  as.matrix()

x_female = data %>%
  filter(sex == "Female") %>%
  select(-country, -year, -Adult, -sex, -Z) %>%
  as.matrix()

id = data %>% select(country, year) %>% distinct() %>% .$country
time = data %>% select(country, year) %>% distinct() %>% .$year

z_male = data %>%
  filter(sex == "Male") %>%
  .$Z %>%
  factor(labels = c("M High", "M High+Adult", "M Low"))

z_female = data %>%
  filter(sex == "Female") %>%
  .$Z %>%
  factor(levels = c("F High", "F Mid", "F Low"))

z_male %>% is.na() %>% sum()
z_female %>% is.na() %>% sum()


#### model fit ####
penal = 1

config = list(
  single_group = FALSE,
  bounds = c(0.01,  5),
  lambda_start = 1,
  n_points = 20,
  n_start = 50,
  n_start_iters = 10,
  n_start_cores = 1,
  epsilon_w = 1e-10,
  beta_sd = sqrt(penal),
  mu_sd = 1,
  sigma_a = 1, sigma_b = 1
)

iters = 500
burn_in = 100
thin = 5
n_cores = 1
n_basis = 10
seed = 1
M = 10
lambda = penal

male_fit2 = fit_mbfseq(
  M = M,
  z = z_male,
  id = id,
  time = time,
  iters = iters,
  burn_in = burn_in,
  thin = thin,
  n_cores = n_cores,
  n_basis = n_basis,
  seed = seed,
  verbose = TRUE,
  config = config,
  lambda = lambda
)


male_fit$metrics
male_fit2$metrics
male_fit2$models$`G=3, M=10`$w_class %>% unique()
male_fit$models$`G=3, M=10`$w_class %>% unique()


male_fit$models$`G=3, M=10`$w_class %>% sort()

pp1 = male_fit$metrics$l_post - male_fit$metrics$l_mean
pp2 = male_fit2$metrics$l_post - male_fit$metrics$l_mean

l1 = male_fit$metrics$l_mean
l2 = male_fit2$metrics$l_mean

pp1
pp2

k1 = (log(27) - 1)*pp1
k2 = (log(27) - 1)*pp2

-2*(l1 - k1)
-2*(l2 - k2)

male_fit %>% plot_seq_cluster(G = 3, M = 10)

male_fit %>%
  posterior_summary(M = 10, G = 3) %>%
  plot_probability

plot(M, -2*male_fit$metrics$l_post + log(27)*(2 * 10 * M + M))


male_fit$metrics$DIC3 %>% plot(x = M, y=.)
male_fit$metrics$pp %>% plot()
2*male_fit$metrics$pp

plot(M, -2*male_fit$metrics$l_mean + (2 * 10 * M + M))


plot(x=M, y=male_fit$metrics$DIC2)
plot(x=M, y=male_fit$metrics$DIC3)

data.frame(
  l1 = male_fit$metrics$l_post,
  l2 = male_fit$metrics$l_mean,
  M = M
) %>%
  gather(l, value, -M) %>%
  ggplot(aes(x=M, y=value, color=l)) +
  geom_line() +
  geom_point()

male_fit %>% plot_seq_cluster(G = 3, M = 6)


z_male %>%
  matrix(nrow = 27, ncol = 60, byrow = TRUE)


female_fit = fit_mbfseq(
  M = M,
  z = z_female,
  id = id,
  time = time,
  iters = iters,
  burn_in = burn_in,
  thin = thin,
  n_cores = n_cores,
  n_basis = n_basis,
  seed = seed,
  verbose = TRUE,
  config = config,
  lambda = lambda
)

plot(M, -2*female_fit$metrics$l_mean + log(27)*(2 * 10 * M + M))


female_fit$metrics$DIC3 %>% plot(x = M, y=.)
female_fit$metrics %>% arrange(DIC3)

male_fit %>% plot_seq_cluster(G = 3, M = 6)
female_fit %>% plot_seq_cluster(G = 3, M = 4)
female_fit %>% plot_seq_cluster(G = 3, M = 5)

data %>%
  filter(sex == "Female") %>%
  mutate(w = female_fit$models$`G=3, M=6`$w_class[country]) %>%
  ggplot(aes(x=`Non-old`, y=Old, color=factor(w))) +
  geom_point()

data %>%
  filter(sex == "Female") %>%
  mutate(w = female_fit$models$`G=3, M=5`$w_class[country]) %>%
  ggplot(aes(x=year, y=`Non-old`, color=factor(w), group=country)) +
  geom_line()



