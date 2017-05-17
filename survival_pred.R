# https://cran.r-project.org/web/packages/survival/survival.pdf
# http://bcb.dfci.harvard.edu/~aedin/courses/Bioconductor/survival.pdf

library(readr)
library(survival)
library(survcomp)

# Load data
X <- read_csv("data/processed/baseline_data_all.csv")

# Prep data - delete extraneous cols, convert outputs to time and event
id <- X$PUBLIC_ID
X$PUBLIC_ID <- NULL

last_obs <-X$LAST_OBSERVED
X$LAST_OBSERVED <- NULL
X$time <- last_obs

cens <- X$CENSORED
X$CENSORED <- NULL
# Surv expects dead = 1/TRUE
dead = cens == 0
X$event <- dead

train <- X$TRAIN_SET == 1
X$TRAIN_SET <- NULL

# Throw out some badly behaved cols
X$DEMOG_AMERICANINDIA <- NULL
X$IC_PATIENTHADANO <- NULL
X$D_IM_reason <- NULL
X$D_TRI_CF_TRISOMIES21 <- NULL

# Train/test split
X_train = X[train,]
X_test = X[!train,]

# Do Cox proportional hazards regression on training set
fit <- coxph(Surv(time, event) ~ ., data = X_train)
summary(fit)

pred_train <- predict(fit, newdata = X_train, type = "risk")
pred_test <- predict(fit, newdata = X_test, type = "risk")

# Concordance index
perf_train <- concordance.index(x = pred_train, surv.time = X_train$time, surv.event = X_train$event, method = "noether", na.rm = TRUE)
print(perf_train[1:5])
perf_test <- concordance.index(x = pred_test, surv.time = X_test$time, surv.event = X_test$event, method = "noether", na.rm = TRUE)
print(perf_test[1:5])
# Overfits

# Time-dependent ROC curve
tau <- 365*1 # 1 yr
perf_train <- tdrocc(x = pred_train, surv.time = X_train$time, surv.event = X_train$event, time = tau,na.rm = TRUE)
perf_tst <- tdrocc(x = pred_test, surv.time = X_test$time, surv.event = X_test$event, time = tau,na.rm = TRUE)
plot(x = 1 - perf_train$spec, y = perf_train$sens, type = "l", lwd = 3, col = "orange",  xlab = "False Positive Rate", ylab = "True Positive Rate", xlim = c(0,1), ylim = c(0,1), main = "Time-dependent ROC curve at 1 yr")
lines(x = 1 - perf_tst$spec, y = perf_tst$sens, type = "l", lwd = 3, col = "blue")
lines(x = c(0,1), y = c(0,1), lty = 3, col = "black")
legend("bottomright", legend=c("Train", "Test"), col=c("orange", "blue"), lty=c(1, 1))
