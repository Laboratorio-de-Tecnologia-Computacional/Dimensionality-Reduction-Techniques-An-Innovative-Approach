## ============================================================
## 0) Full reproducibility + utilities (journal-ready)
## ============================================================

# Pin threads (determinism across BLAS/OMP implementations)
Sys.setenv(OMP_NUM_THREADS = "1",
           OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1",
           NUMEXPR_NUM_THREADS = "1")

# Deterministic seed (stable across machines)
set.seed(123, kind = "Mersenne-Twister", normal.kind = "Inversion")

`%||%` <- function(x, y) if (is.null(x) || length(x)==0) y else x

# --- Helper functions ---
# Extracts 2D coordinates from dimRedResult or from a raw matrix
get_coords2d <- function(obj) {
  if (inherits(obj, "dimRedResult")) as.matrix(obj@data@data)[, 1:2, drop = FALSE]
  else as.matrix(obj)[, 1:2, drop = FALSE]
}

# Attempts to infer a class/label column for coloring the scatter plot
infer_label_vec <- function(df_raw) {
  cand <- c("Diagnosis","Part","class","Class","Type","target","Label","Classe","Grupo","Group")
  nm <- cand[cand %in% names(df_raw)]
  if (length(nm)) return(as.factor(df_raw[[nm[1]]]))
  factor(rep("All", nrow(df_raw)))
}

# k-NN indices (excludes the point itself)
get_knn_idx <- function(X, k) {
  idx <- RANN::nn2(X, query = X, k = k + 1)$nn.idx
  idx[, -1, drop = FALSE]
}

# R_NX curve (average neighbor overlap for k = 1..kmax)
rn_curve <- function(X_high, X_low, kmax) {
  N <- nrow(X_high)
  nn_high_all <- get_knn_idx(X_high, kmax)
  nn_low_all  <- get_knn_idx(X_low,  kmax)
  vapply(1:kmax, function(k) {
    nh_k <- nn_high_all[, 1:k, drop = FALSE]
    nl_k <- nn_low_all[,  1:k, drop = FALSE]
    mean(vapply(1:N, function(i) length(intersect(nh_k[i,], nl_k[i,]))/k, numeric(1)))
  }, numeric(1))
}

# Normalized trapezoidal AUC for a sequence y[1..kmax]
auc_trapz_norm <- function(y) {
  if (length(y) < 2) return(NA_real_)
  auc <- sum((y[-1] + y[-length(y)]) / 2)
  auc / (length(y) - 1)
}

# Trustworthiness/Continuity fallback (pure R), plus Precision/Recall proxies
trust_continuity <- function(X_high, X_low, k) {
  n <- nrow(X_high)
  k <- min(k, n - 1L)
  
  # Compute pairwise distance matrices
  D_high <- as.matrix(dist(X_high))
  D_low  <- as.matrix(dist(X_low))
  
  # Ranking of neighbors (column i = ranking of other points relative to point i)
  ord_high <- apply(D_high, 1, order)
  ord_low  <- apply(D_low,  1, order)
  
  # k-nearest neighbors (excluding the point itself at position 1)
  N_high <- ord_high[2:(k + 1), , drop = FALSE]  # k × n matrix
  N_low  <- ord_low[ 2:(k + 1), , drop = FALSE]  # k × n matrix
  
  trust_sum <- 0
  cont_sum  <- 0
  prec_i    <- numeric(n)
  
  for (i in 1:n) {
    # Intrusions: neighbors present in the low-dimensional space but not in high-dimensional space
    U <- setdiff(N_low[, i], N_high[, i])
    if (length(U) > 0) {
      ranks_high <- match(U, ord_high[, i])   # ranks 1..n
      trust_sum  <- trust_sum + sum(pmax(0, ranks_high - k))
    }
    
    # Misses: neighbors present in the high-dimensional space but not in low-dimensional space
    V <- setdiff(N_high[, i], N_low[, i])
    if (length(V) > 0) {
      ranks_low <- match(V, ord_low[, i])
      cont_sum  <- cont_sum + sum(pmax(0, ranks_low - k))
    }
    
    # Precision proxy: intersection between high-dim and low-dim neighborhoods
    inter <- intersect(N_low[, i], N_high[, i])
    prec_i[i] <- length(inter) / k  # normalized (0–1)
  }
  
  denom <- n * k * (2 * n - 3 * k - 1)
  
  # If the denominator is invalid (k too large), reduce k and recompute
  if (denom <= 0) {
    k2 <- max(1L, floor((2 * n - 2) / 3))
    if (k2 < k) {
      return(trust_continuity(X_high, X_low, k2))
    }
  }
  
  # Trustworthiness and Continuity formulas
  T <- 1 - 2 * trust_sum / denom
  C <- 1 - 2 * cont_sum  / denom
  
  # Clamp values into the [0, 1] range
  T <- max(min(T, 1), 0)
  C <- max(min(C, 1), 0)
  
  list(
    Trustworthiness = T,
    Continuity      = C,
    Precision       = mean(prec_i),
    Recall          = mean(prec_i)
  )
}


## ============================================================
## 1) DATA HUB — unified loaders and preprocessing for 7 datasets
## ============================================================

# -- df1: Breast Cancer Wisconsin (Diagnostic) --
load_df1 <- function() {
  url1 <- "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data"
  df1 <- read.csv(url1, header = FALSE)
  colnames(df1) <- c(
    "ID", "Diagnosis",
    "radius_mean","texture_mean","perimeter_mean","area_mean","smoothness_mean",
    "compactness_mean","concavity_mean","concave_points_mean","symmetry_mean","fractal_dimension_mean",
    "radius_se","texture_se","perimeter_se","area_se","smoothness_se",
    "compactness_se","concavity_se","concave_points_se","symmetry_se","fractal_dimension_se",
    "radius_worst","texture_worst","perimeter_worst","area_worst","smoothness_worst",
    "compactness_worst","concavity_worst","concave_points_worst","symmetry_worst","fractal_dimension_worst"
  )
  df1 <- df1[, -1]  # drop ID
  df1
}

# -- df2: Chemical Composition of Ceramic --
load_df2 <- function(path_csv =
                       "C:/Users/Usuario/Desktop/RD - Trabalho/Artigo Completo/Datasets/2-Chemical Composition of Ceramic/Chemical Composion of Ceramic.csv") { 
  df2 <- read.csv(path_csv, check.names = FALSE)
  if (names(df2)[1] %in% c("ID","Id","id","#","index")) df2 <- df2[, -1]
  cols_wt <- c("Na2O","MgO","Al2O3","SiO2","K2O","CaO","TiO2","Fe2O3")
  if (all(cols_wt %in% names(df2))) {
    df2 <- df2 %>%
      dplyr::mutate(across(all_of(cols_wt), ~ .x * 10000)) %>%  # %wt → ppm
      dplyr::mutate(across(where(is.numeric), abs))
  } else {
    df2 <- df2 %>% dplyr::mutate(across(where(is.numeric), abs))
  }
  df2
}

# -- df3: Dry Bean (Excel) --
load_df3 <- function(path_xlsx =
                       "C:/Users/Usuario/Desktop/RD - Trabalho/Artigo Completo/Datasets/4-Dry Bean/Dry_Bean_Dataset.xlsx") {
  df3 <- readxl::read_excel(path_xlsx)
  if ("Class" %in% names(df3)) df3 <- df3 %>% rename(class = Class)
  df3$class <- as.factor(df3$class)
  as.data.frame(df3)
}

# -- df4: Glass Identification (UCI) --
load_df4 <- function() {
  url4 <- "https://archive.ics.uci.edu/ml/machine-learning-databases/glass/glass.data"
  df4 <- read.csv(url4, header = FALSE)
  colnames(df4) <- c("Id","RI","Na","Mg","Al","Si","K","Ca","Ba","Fe","Type")
  df4 <- df4[, -1]
  df4$Type <- factor(df4$Type,
                     levels = c(1,2,3,5,6,7),
                     labels = c("bw_float","bw_non_float","vehicle_windows","containers","tableware","headlamps"))
  df4
}

# -- df5: Hepatitis (UCI) --
load_df5 <- function() {
  url5 <- "https://archive.ics.uci.edu/ml/machine-learning-databases/hepatitis/hepatitis.data"
  df5 <- read.csv(url5, header = FALSE, na.strings = "?")
  colnames(df5) <- c("Class","AGE","SEX","STEROID","ANTIVIRALS","FATIGUE","MALAISE","ANOREXIA",
                     "LIVER_BIG","LIVER_FIRM","SPLEEN_PALPABLE","SPIDERS","ASCITES","VARICES",
                     "BILIRUBIN","ALK_PHOSPHATE","SGOT","ALBUMIN","PROTIME","HISTOLOGY")
  df5 <- df5 %>%
    mutate(
      Class = factor(Class, levels = c(1,2), labels = c("Die","Live")),
      SEX = factor(SEX, levels = c(1,2), labels = c("male","female")),
      across(c(STEROID,ANTIVIRALS,FATIGUE,MALAISE,ANOREXIA,LIVER_BIG,LIVER_FIRM,
               SPLEEN_PALPABLE,SPIDERS,ASCITES,VARICES,HISTOLOGY),
             ~ factor(., levels=c(1,2), labels=c("no","yes"))),
      across(c(BILIRUBIN,ALK_PHOSPHATE,SGOT,ALBUMIN,PROTIME,AGE), as.numeric)
    )
  num_cols <- names(df5)[sapply(df5, is.numeric)]
  fac_cols <- names(df5)[sapply(df5, is.factor)]
  if (length(num_cols))
    df5[num_cols] <- lapply(df5[num_cols], function(x){ x[is.na(x)] <- median(x, na.rm=TRUE); x })
  if (length(fac_cols))
    df5[fac_cols] <- lapply(df5[fac_cols], function(x){
      m <- names(which.max(table(x))); y <- as.character(x); y[is.na(y)] <- m; factor(y)
    })
  df5
}

# -- df6: Iris (UCI) --
load_df6 <- function() {
  url6 <- "https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data"
  df6 <- read.csv(url6, header = FALSE, na.strings = c("", "?"))
  colnames(df6) <- c("SepalLength","SepalWidth","PetalLength","PetalWidth","Class")
  df6 <- df6 %>% filter(!if_all(everything(), ~ is.na(.)))
  df6$Class <- factor(df6$Class,
                      levels = c("Iris-setosa","Iris-versicolor","Iris-virginica"),
                      labels = c("setosa","versicolor","virginica"))
  df6
}

# -- df7: Wine (UCI) --
load_df7 <- function() {
  url7 <- "https://archive.ics.uci.edu/ml/machine-learning-databases/wine/wine.data"
  df7 <- read.csv(url7, header = FALSE, na.strings = c("", "?"))
  colnames(df7) <- c("Class","Alcohol","Malic_Acid","Ash","Alcalinity_of_Ash","Magnesium",
                     "Total_Phenols","Flavanoids","Nonflavanoid_Phenols","Proanthocyanins",
                     "Color_Intensity","Hue","OD280_OD315","Proline")
  df7 <- df7 %>% filter(!if_all(everything(), ~ is.na(.)))
  df7
}

# --- Select the active dataset (single switch for the entire pipeline) ---
ACTIVE_DATASET <- "df3"  # choices: "df1","df2","df3","df4","df5","df6","df7"

# Load all (explicit for review reproducibility)
df1 <- try(load_df1(), silent = TRUE)
df2 <- try(load_df2(), silent = TRUE)
df3 <- try(load_df3(), silent = TRUE)
df4 <- try(load_df4(), silent = TRUE)
df5 <- try(load_df5(), silent = TRUE)
df6 <- try(load_df6(), silent = TRUE)
df7 <- try(load_df7(), silent = TRUE)


# Bind the selected dataset to generic names used downstream
df <- switch(ACTIVE_DATASET,
             df1 = { if (inherits(df1,"try-error")) stop(attr(df1,"condition")$message); df1 },
             df2 = { if (inherits(df2,"try-error")) stop(attr(df2,"condition")$message); df2 },
             df3 = { if (inherits(df3,"try-error")) stop(attr(df3,"condition")$message); df3 },
             df4 = { if (inherits(df4,"try-error")) stop(attr(df4,"condition")$message); df4 },
             df5 = { if (inherits(df5,"try-error")) stop(attr(df5,"condition")$message); df5 },
             df6 = { if (inherits(df6,"try-error")) stop(attr(df6,"condition")$message); df6 },
             df7 = { if (inherits(df7,"try-error")) stop(attr(df7,"condition")$message); df7 },
             stop("ACTIVE_DATASET must be one of: df1..df7")
)

# Standardized objects for the DR pipeline
label_vec <- infer_label_vec(df)
df_num    <- df[, sapply(df, is.numeric), drop = FALSE]
df_scaled <- scale(df_num)
X_high    <- as.matrix(df_scaled)
N         <- nrow(X_high)

## ============================================================
## 2) Dimensionality reduction — methods, seeds, deterministic execution
## ============================================================

global_methods <- c("DRR","Isomap","PCA","PCA_L1","kPCA","FastICA","MDS","nMDS","DrL","FruchtermanReingold","KamadaKawai")
local_methods  <- c("HLLE","DiffusionMaps","tSNE")
methods <- unique(c(global_methods, local_methods))

method_pkgs <- list(
  DiffusionMaps="diffusionMap", DRR="DRR", FastICA="fastICA",
  KamadaKawai="igraph", DrL="igraph", FruchtermanReingold="igraph",
  HLLE=c("RANN","RSpectra","Matrix"), Isomap="vegan", kPCA="kernlab",
  PCA_L1="pcaL1", MDS="stats", nMDS="MASS", PCA="stats", tSNE="Rtsne"
)

pkg_ok <- function(pkgs) all(sapply(pkgs, function(p) if (p=="stats") TRUE else requireNamespace(p, quietly=TRUE)))
methods_ok <- methods[sapply(methods, function(m) pkg_ok(method_pkgs[[m]]))]

data_dimred <- dimRedData(as.data.frame(df_scaled))

# Stable per-method seed (simple hash without external deps)
seed_for_method <- function(m) sum(utf8ToInt(m)) + 2024L

run_dimred <- function(m) {
  set.seed(seed_for_method(m), kind="Mersenne-Twister", normal.kind="Inversion")
  if (m == "tSNE") {
    per_max <- max(5, floor((N - 1) / 3) - 1)
    per_use <- max(5, min(30, per_max))
    res <- tryCatch(embed(data_dimred, "tSNE", ndim=2,
                          perplexity=per_use, theta=.5, max_iter=1000, pca=FALSE),
                    error=function(e) NULL)
    return(res)
  } else {
    return(tryCatch(embed(data_dimred, m, ndim=2), error=function(e) NULL))
  }
}

dimred_results <- purrr::map(methods_ok, run_dimred) %>% rlang::set_names(methods_ok)
dimred_results <- dimred_results[!vapply(dimred_results, is.null, TRUE)]

## ============================================================
## 3) Metrics (13+)
## ============================================================

# 3.1 Q_local / Q_global (dimRed internal)
qs <- purrr::map_dfr(names(dimred_results), function(m){
  proj <- dimred_results[[m]]
  ql <- tryCatch(dimRed::quality(proj, "Q_local"),  error=function(e) NA_real_)
  qg <- tryCatch(dimRed::quality(proj, "Q_global"), error=function(e) NA_real_)
  tibble::tibble(Method=m, Q_local=ql, Q_global=qg,
                 Group=dplyr::case_when(m %in% global_methods ~ "Global",
                                        m %in% local_methods  ~ "Local",
                                        TRUE ~ "Other"))
})

# 3.2 R_NX curve + normalized AUC
kmax <- min(30, N-1)
curves_list <- list(); auc_tbl <- tibble::tibble()
for (m in names(dimred_results)) {
  coords <- tryCatch(get_coords2d(dimred_results[[m]]), error=function(e) NULL)
  if (is.null(coords) || nrow(coords)!=N || ncol(coords)<2) next
  Rk <- rn_curve(X_high, coords[,1:2, drop=FALSE], kmax)
  curves_list[[m]] <- Rk
  auc_tbl <- dplyr::bind_rows(auc_tbl, tibble::tibble(Method=m, AUC=auc_trapz_norm(Rk)))
}
auc_tbl <- dplyr::arrange(auc_tbl, dplyr::desc(AUC))  # AUC_R_NX downstream

# 3.3 DRquality + fallbacks
Data <- X_high
Cls_int <- if (nlevels(label_vec) > 1) as.integer(label_vec) else rep(1L, N)

extract_measures <- function(m) {
  set.seed(seed_for_method(m))
  proj <- tryCatch(get_coords2d(dimred_results[[m]]), error=function(e) NULL)
  if (is.null(proj) || nrow(proj)!=nrow(Data)) return(NULL)
  k_use <- min(15, nrow(Data) - 1)
  
  td <- tryCatch(DRquality::MeasureTandD(Data, pData=proj, NeighborhoodSize=k_use), error=function(e) NULL)
  pr <- tryCatch(DRquality::PrecisionAndRecall(Data, pData=proj, NeighborhoodSize=k_use), error=function(e) NULL)
  
  if (is.null(td) || is.null(pr)) {
    tc <- trust_continuity(Data, proj, k_use)
    trust_mean <- as.numeric(tc$Trustworthiness)
    cont_mean  <- as.numeric(tc$Continuity)
    disc_mean  <- 1 - cont_mean
    prec_mean  <- as.numeric(tc$Precision)
    rec_mean   <- as.numeric(tc$Recall)
  } else {
    get_mean <- function(x, name, idx) {
      if (is.null(x)) return(NA_real_)
      if (!is.null(dim(x))) { xd <- as.data.frame(x); if (name %in% names(xd)) return(mean(xd[[name]], na.rm=TRUE)); if (ncol(xd)>=idx) return(mean(xd[[idx]], na.rm=TRUE)) }
      if (is.list(x)) { if (!is.null(x[[name]])) return(mean(x[[name]], na.rm=TRUE)); if (!is.null(x[[idx]])) return(mean(x[[idx]], na.rm=TRUE)) }
      if (is.atomic(x)) { if (!is.null(names(x)) && name %in% names(x)) return(as.numeric(x[[name]])); if (length(x)>=idx) return(as.numeric(x[[idx]])) }
      NA_real_
    }
    trust_mean <- get_mean(td, "Trustworthiness", 1)
    cont_mean  <- get_mean(td, "Continuity",      2)
    disc_mean  <- if (is.finite(cont_mean)) 1 - cont_mean else get_mean(td, "Discontinuity", 2)
    prec_mean  <- get_mean(pr, "Precision", 1)
    rec_mean   <- get_mean(pr, "Recall",    2)
  }
  
  kt  <- tryCatch(DRquality::KendallsTau(as.matrix(dist(Data)), as.matrix(dist(proj))), error=function(e) NA_real_)
  ce  <- tryCatch(DRquality::ClassificationError(as.matrix(dist(proj)), Cls_int, k=5), error=function(e) list(Accuracy=NA))
  gce <- tryCatch(DRquality::GabrielClassificationError(Data, proj, Cls_int, PlotIt=FALSE), error=function(e) list(GCE=NA))
  cm  <- tryCatch(DRquality::Cmeasure(Data, proj), error=function(e) NA_real_)
  acc_val <- if (is.list(ce)  && !is.null(ce$Accuracy)) ce$Accuracy else as.numeric(ce)[1]
  gce_val <- if (is.list(gce) && !is.null(gce$GCE))      gce$GCE      else as.numeric(gce)[1]
  
  tibble::tibble(
    Method=m, Trustworthiness=trust_mean, Continuity=cont_mean, Discontinuity=disc_mean,
    Precision=prec_mean,     Recall=rec_mean, KendallTau=as.numeric(kt)[1],
    Accuracy_kNN=acc_val,    GCE=gce_val,     C_Measure=as.numeric(cm)[1]
  )
}
drq_results <- purrr::map_dfr(names(dimred_results), extract_measures)

# 3.4 Cophenetic correlation (distance preservation)
calc_cophenetic_all <- function(dimred_results, X_high) {
  purrr::map_dfr(names(dimred_results), function(m){
    coords <- tryCatch(get_coords2d(dimred_results[[m]]), error=function(e) NULL)
    if (is.null(coords) || nrow(coords)!=nrow(X_high)) return(tibble::tibble(Method=m, Cophenetic=NA_real_))
    r <- suppressWarnings(cor(as.vector(dist(X_high)), as.vector(dist(coords)), method="pearson"))
    tibble::tibble(Method=m, Cophenetic=as.numeric(r))
  })
}
cophen_tbl <- calc_cophenetic_all(dimred_results, X_high)

# 3.5 Reconstruction RMSE (only for methods with inverse map)
methods_with_recon <- c("PCA","FastICA","kPCA","DRR")
rmse_tbl <- purrr::map_dfr(methods_with_recon, function(m){
  if (is.null(dimred_results[[m]])) return(tibble::tibble(Method=m, RMSE_Reconstruction=NA_real_))
  set.seed(seed_for_method(m))
  val <- tryCatch(dimRed::reconstruction_rmse(dimred_results[[m]]), error=function(e) NA_real_)
  tibble::tibble(Method=m, RMSE_Reconstruction=as.numeric(val))
})

## ============================================================
## 4) Unified metrics table
## ============================================================
qs_tbl   <- dplyr::select(qs, Method, Q_local, Q_global)
auc_tbl2 <- dplyr::rename(auc_tbl, AUC_R_NX = AUC)

# Merge all metrics into a single table (initial version)
all_metrics <- drq_results %>%
  dplyr::full_join(qs_tbl,     by = "Method") %>%
  dplyr::full_join(auc_tbl2,   by = "Method") %>%
  dplyr::full_join(cophen_tbl, by = "Method") %>%
  dplyr::full_join(rmse_tbl,   by = "Method")

## 4.1 Recompute Trustworthiness / Continuity / Precision / Recall for ALL methods
k_use <- min(15, nrow(Data) - 1L)

trust_tbl <- purrr::map_dfr(names(dimred_results), function(m) {
  proj <- tryCatch(get_coords2d(dimred_results[[m]]), error = function(e) NULL)
  
  # Handle missing projections
  if (is.null(proj) || nrow(proj) != nrow(Data)) {
    return(tibble::tibble(
      Method          = m,
      Trustworthiness = NA_real_,
      Continuity      = NA_real_,
      Precision       = NA_real_,
      Recall          = NA_real_
    ))
  }
  
  # Compute trustworthiness metrics using the custom function
  tc <- trust_continuity(Data, proj, k_use)
  
  tibble::tibble(
    Method          = m,
    Trustworthiness = as.numeric(tc$Trustworthiness),
    Continuity      = as.numeric(tc$Continuity),
    Precision       = as.numeric(tc$Precision),
    Recall          = as.numeric(tc$Recall)
  )
})

# Remove old columns (if present) and append the recomputed metrics
all_metrics <- all_metrics %>%
  dplyr::select(-dplyr::any_of(c("Trustworthiness", "Continuity",
                                 "Precision", "Recall"))) %>%
  dplyr::left_join(trust_tbl, by = "Method") %>%
  
  # Reorder columns for cleaner output
  dplyr::relocate(Method, Q_local, Q_global, AUC_R_NX, Cophenetic, RMSE_Reconstruction,
                  Trustworthiness, Continuity, Discontinuity, Precision, Recall,
                  KendallTau, Accuracy_kNN, GCE, C_Measure)

# Inspect updated metrics
all_metrics %>%
  dplyr::select(Method, Trustworthiness, Continuity, Precision, Recall)

## 4.2 Correlation matrix and corrplot
is_num   <- vapply(all_metrics, is.numeric, logical(1))
mat_num  <- all_metrics[ , is_num, drop = FALSE]

# Pairwise correlations across all available values
cor_mat <- stats::cor(mat_num, use = "pairwise.complete.obs")

# Replace NA / NaN / Inf with 0 to avoid missing tiles in corrplot
cor_mat[!is.finite(cor_mat)] <- 0

# Generate the correlogram
corrplot::corrplot(cor_mat, method = "color", tl.cex = 0.6)

## ============================================================
## 5) Objective-based ranking (baseline + ±10% robustness)
## ============================================================
larger_better  <- c("Q_local","Q_global","AUC_R_NX","Cophenetic",
                    "Trustworthiness","Precision","Recall","KendallTau","Accuracy_kNN")
smaller_better <- c("RMSE_Reconstruction","Discontinuity","GCE","C_Measure")

norm_01 <- function(x, larger=TRUE){
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  rg <- range(x, na.rm=TRUE); s <- if (diff(rg)==0) rep(.5, length(x)) else (x - rg[1]) / diff(rg)
  if (!larger) s <- 1 - s; s
}
rank_dir <- function(x, larger=TRUE){
  r <- if (larger) rank(-x, na.last="keep", ties.method="min") else rank(x, na.last="keep", ties.method="min")
  r[is.na(r)] <- max(r, na.rm=TRUE) + 1L
  as.integer(r)
}

get_scenario_weights <- function(scenario = c("clustering","reconstruction","visualization","ml_pipeline")) {
  scenario <- match.arg(scenario)
  switch(scenario,
         clustering     = c(Local=.50, Global=.20, Geom=.15, Class=.15),
         reconstruction = c(Local=.30, Global=.30, Geom=.40, Class=.00),
         visualization  = c(Local=.40, Global=.20, Geom=.00, Class=.40),
         ml_pipeline    = c(Local=.40, Global=.20, Geom=.00, Class=.40)
  )
}

perturb_weights <- function(w, delta = 0.10) {
  bases <- c(Local = w["Local"], 
             Global = w["Global"], 
             Geom = w["Geom"], 
             Class = w["Class"])
  mults <- expand.grid(
    Local  = bases["Local"]  * c(-delta, +delta),
    Global = bases["Global"] * c(-delta, +delta),
    Geom   = bases["Geom"]   * c(-delta, +delta),
    Class  = bases["Class"]  * c(-delta, +delta),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )
  ws <- apply(mults, 1, function(row){
    v <- bases + row
    v <- pmax(v, 0) 
    v / sum(v)      
  })
  
  lapply(seq_len(ncol(ws)), function(i) {
    x <- ws[, i]
    names(x) <- names(bases)
    x
  })
}
  
compute_indices <- function(df_norm) {
  need <- function(nms) rowMeans(cbind(df_norm[, intersect(nms, names(df_norm)), drop=FALSE]), na.rm=TRUE)
  df_idx <- df_norm
  df_idx$Index_Local  <- need(c("Q_local","Trustworthiness","Precision","Recall"))
  df_idx$Index_Global <- need(c("Q_global","AUC_R_NX","Cophenetic","KendallTau"))
  df_idx$Index_Geom   <- need(c("RMSE_Reconstruction","Discontinuity","GCE","C_Measure"))
  df_idx$Index_Class  <- need(c("Accuracy_kNN","Precision","Recall"))
  df_idx
}

normalize_metrics_table <- function(df) {
  out <- df
  mets <- setdiff(names(df), "Method")
  for (m in mets) {
    if (m %in% larger_better) out[[m]] <- norm_01(df[[m]], TRUE)
    else if (m %in% smaller_better) out[[m]] <- norm_01(df[[m]], FALSE)
  }
  out
}

best_metric_and_champion <- function(all_metrics) {
  metrics_all <- intersect(c(larger_better, smaller_better), names(all_metrics))
  mat <- all_metrics %>% dplyr::select(Method, dplyr::all_of(metrics_all))
  scores_long <- mat %>%
    tidyr::pivot_longer(-Method, names_to="Metric", values_to="Value") %>%
    group_by(Metric) %>%
    mutate(Larger = Metric %in% larger_better,
           Rank   = rank_dir(Value, unique(Larger))) %>%
    ungroup()
  rank_wide <- scores_long %>%
    dplyr::select(Method, Metric, Rank) %>%
    tidyr::pivot_wider(names_from = Metric, values_from = Rank)
  rownames(rank_wide) <- rank_wide$Method
  rank_wide <- rank_wide[, setdiff(names(rank_wide), "Method"), drop=FALSE]
  
  best_metric <- NA_character_
  if (ncol(rank_wide) >= 2) {
    cor_rank <- suppressWarnings(cor(rank_wide, use="pairwise.complete.obs", method="spearman"))
    cor_abs  <- abs(cor_rank); diag(cor_abs) <- NA
    cons_metric <- tibble::tibble(
      Metric = colnames(cor_abs),
      Consistency = apply(cor_abs, 2, function(col) mean(col, na.rm=TRUE))
    ) %>% arrange(desc(Consistency))
    best_metric <- cons_metric$Metric[1]
  }
  
  champion <- NA_character_
  if (!is.na(best_metric)) {
    colv <- all_metrics[[best_metric]]
    champion <- if (best_metric %in% larger_better)
      all_metrics$Method[which.max(colv)] else all_metrics$Method[which.min(colv)]
  }
  list(best_metric = best_metric, champion = champion)
}

rank_by_objective <- function(all_metrics, scenario) {
  df_norm <- normalize_metrics_table(all_metrics)
  df_idx  <- compute_indices(dplyr::rename(df_norm, algorithm=Method))
  w <- get_scenario_weights(scenario)
  
  df_rank <- df_idx %>%
    mutate(Final_Score = w["Local"]*Index_Local + w["Global"]*Index_Global +
             w["Geom"] *Index_Geom  + w["Class"] *Index_Class) %>%
    arrange(desc(Final_Score))
  df_rank$Rank <- seq_len(nrow(df_rank))
  
  w_list <- perturb_weights(w, delta=0.10)
  ranks_pert <- lapply(w_list, function(wp){
    sc <- wp["Local"]*df_idx$Index_Local + wp["Global"]*df_idx$Index_Global +
      wp["Geom"] *df_idx$Index_Geom  + wp["Class"] *df_idx$Index_Class
    ord <- order(-sc, df_idx$algorithm)
    data.frame(algorithm = df_idx$algorithm[ord],
               Score = sc[ord],
               Rank = seq_along(sc),
               stringsAsFactors = FALSE)
  })
  top1 <- sapply(ranks_pert, function(x) x$algorithm[1])
  top1_freq <- table(top1)
  top1_freq <- tibble::tibble(algorithm=names(top1_freq),
                              Top1_Freq = as.integer(top1_freq),
                              Top1_Freq_Pct = as.numeric(top1_freq)/length(ranks_pert))
  score_mean <- Reduce(function(a,b) merge(a,b,by="algorithm",all=TRUE),
                       lapply(ranks_pert, function(x) x[,c("algorithm","Score")]))
  score_mean$Score_Mean_Pert <- rowMeans(score_mean[,-1], na.rm=TRUE)
  score_mean <- score_mean[,c("algorithm","Score_Mean_Pert")]
  
  mmc <- best_metric_and_champion(all_metrics)
  
  list(
    ranking_baseline = df_rank,
    robustness_top1  = dplyr::arrange(dplyr::left_join(top1_freq, score_mean, by="algorithm"),
                                      dplyr::desc(Top1_Freq), dplyr::desc(Score_Mean_Pert)),
    best_metric      = mmc$best_metric,
    champion_metric  = mmc$champion
  )
}

## ============================================================
## 6) Run objectives + plots of the champion under the best metric
## ============================================================
metric_names_pretty <- c(
  Q_local="Q_local", Q_global="Q_global", AUC_R_NX="AUC(R_NX)",
  Cophenetic="Cophenetic correlation", Trustworthiness="Trustworthiness",
  Precision="Precision", Recall="Recall", KendallTau="Kendall's Tau",
  Accuracy_kNN="kNN Accuracy", RMSE_Reconstruction="Reconstruction RMSE",
  Discontinuity="Discontinuity", GCE="Gabriel CE", C_Measure="C-measure"
)
                              
save_plot <- function(gg, filename, width=9, height=8, dpi=300) {
  ggsave(filename, gg, width=width, height=height, dpi=dpi)
}

plot_champion_best_metric <- function(best_metric, best_method, objective_name) {
  stopifnot(!is.null(best_metric), !is.null(best_method), nzchar(best_metric), nzchar(best_method))
  obj <- dimred_results[[best_method]]
  stopifnot(!is.null(obj))
  coords <- get_coords2d(obj); stopifnot(is.matrix(coords), ncol(coords)>=2)
  
  val <- tryCatch(dplyr::filter(all_metrics, .data$`Method`==best_method)[[best_metric]], error=function(e) NA_real_)
  val_txt <- if (length(val) && is.finite(val[1])) paste0(" — value = ", signif(val[1],5)) else ""
  df_plot <- data.frame(x=coords[,1], y=coords[,2], Label=label_vec)
  title_metric <- if (!is.null(metric_names_pretty[[best_metric]])) metric_names_pretty[[best_metric]] else best_metric
  
  gg <- ggplot(df_plot, aes(x=x, y=y, color=Label)) +
    geom_point(alpha=1, size=4) +
    theme_minimal(base_size=19.5) +
    labs(
      title=paste0("Objective: ", objective_name, " — Champion of the best metric: ", title_metric),
      subtitle=paste0(best_method, val_txt),
      x="Dim 1", y="Dim 2", color="Class"
    ) +
    theme(plot.title=element_text(face="bold"))
  print(gg)
  #invisible(gg)
                  
 fname <- paste0("plot_", objective_name, "_best_metric_", best_metric, "_", best_method, ".png")
 save_plot(gg, fname)
}

objectives <- c("clustering","reconstruction","visualization","ml_pipeline")

obj_results <- vector("list", length(objectives))
names(obj_results) <- objectives

objectives_scores_long <- vector("list", length(objectives))
names(objectives_scores_long) <- objectives

winners_by_metric <- NULL
winners_by_objective <- NULL


# Winners per metric (extreme on the raw metric)
for (met in c(larger_better, smaller_better)) {
  if (!met %in% names(all_metrics)) next
  colv <- all_metrics[[met]]
  best_method <- if (met %in% larger_better) {
    all_metrics$Method[which.max(colv)]
  } else {
    all_metrics$Method[which.min(colv)]
  }
  winners_by_metric <- dplyr::bind_rows(
    winners_by_metric,
    tibble::tibble(Metric = met, Best_Method = best_method,
                   Value = as.numeric(all_metrics[all_metrics$Method==best_method, met]))
  )
}

plot_best_method_by_objective <- function(best_method, objective_name) {
  stopifnot(!is.null(best_method), nzchar(best_method))
  obj <- dimred_results[[best_method]]
  stopifnot(!is.null(obj))
  
  coords <- get_coords2d(obj)
  stopifnot(is.matrix(coords), ncol(coords)>=2)
  
  df_plot <- data.frame(x = coords[,1], y = coords[,2], Label = label_vec)
  
  gg <- ggplot(df_plot, aes(x=x, y=y, color=Label)) +
    geom_point(alpha=1, size=4) +
    theme_minimal(base_size=19.5) +
    labs(
      title = paste0("Objective: ", objective_name,
                     " — Champion by Final Score (Baseline Ranking)"),
      subtitle = best_method,
      x="Dim 1", y="Dim 2", color="Class"
    ) +
    theme(plot.title = element_text(face="bold"))
  
  print(gg)
 fname <- paste0("plot_", objective_name, "_champion_baseline_", best_method, ".png")
 save_plot(gg, fname)
}

objectives_scores_long <- list()
winners_by_objective <- tibble::tibble()

for (obj in objectives) {
  cat("\n================= OBJECTIVE:", obj, "=================\n")
  res <- rank_by_objective(all_metrics, obj)
  obj_results[[obj]] <- res
  
  print(head(res$ranking_baseline, 10))
  best_method_baseline <- res$ranking_baseline$algorithm[1]
  plot_best_method_by_objective(best_method_baseline, obj)
  cat("Best metric (consistency):", res$best_metric, "\n")
  cat("Champion under that metric:", res$champion_metric, "\n")
  
  if (!is.na(res$best_metric) && !is.na(res$champion_metric)) {
    plot_champion_best_metric(res$best_metric, res$champion_metric, obj)
  }
  
  # Baseline ranking bar plot
  gg <- ggplot(res$ranking_baseline, aes(x=reorder(algorithm, Final_Score), y=Final_Score)) +
    geom_col(alpha=1,  fill="black") + coord_flip() +
    geom_text(aes(label=sprintf("%.3f", Final_Score)), hjust=-.1, size=8) +
    labs(title=paste0("DR ranking — Objective: ", obj), x="Algorithm", y="Score (0–1)") +
    theme_minimal(base_size=22) +
    ylim(0, max(res$ranking_baseline$Final_Score, na.rm=TRUE)*1.1)
  print(gg)
  
  # Long table for Excel (baseline + robustness summary)
  
  objectives_scores_long[[obj]] <- res$ranking_baseline %>%
    dplyr::select(algorithm,
                  Index_Local, Index_Global, Index_Geom, Index_Class,
                  Final_Score) %>%
    dplyr::mutate(Objective = obj)
  
  winners_by_objective <- dplyr::bind_rows(
    winners_by_objective,
    tibble::tibble(
      Objective              = obj,
      Best_Metric            = res$best_metric,
      Champion_of_Best_Metric = res$champion_metric,
      Top1_Baseline          = res$ranking_baseline$algorithm[1]
    ))
}

objectives_scores_long <- dplyr::bind_rows(objectives_scores_long)

sapply(obj_results, function(x) x$best_metric)
## ============================================================
## 7) Export Excel with four sheets
## ============================================================
writexl::write_xlsx(
  list(
    metrics_all           = all_metrics,
    objectives_scores     = objectives_scores_long,
    winners_per_metric    = winners_by_metric,
    winners_per_objective = winners_by_objective
  ),
  path = "DR_complete_results.xlsx"
)
message("Excel saved: DR_complete_results.xlsx (4 sheets).")



