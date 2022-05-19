



cif.se <- function(time, ci, n.risk, n.event, kmsurv, kmsummand){
  # from https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Cumulative_Incidence.pdf
  # also here: https://onlinelibrary.wiley.com/doi/full/10.1002/bimj.200900039?saml_referrer
  timeindex <- seq_along(time)
  lapply(timeindex, function(i) {
    cii <- ci[1:i]
    bi <- ((n.risk - n.event)/n.risk)[1:i]
    di <- (n.event/(n.risk^2))[1:i]
    lagkmi <- lag(kmsurv,1,1)[1:i]
    kmsummandi <- kmsummand[1:i]
    vt <-
      sum(((cii[i] - cii)^2) * kmsummandi) +
      sum((lagkmi^2) * bi * di) +
      -2* sum((cii[i] - cii) * lagkmi * di)

    sqrt(vt)
  })
}



testdata = data_matched %>%
  mutate(
    statuskm = fct_recode(status, "censored"="noncoviddeath", "censored"="coviddeath"),
    statusci = status
  )


testfitkm <- survfit(Surv(tte_outcome, ind_outcome) ~ 1, data = testdata, conf.type="plain")
testfitstatuskm <- survfit(Surv(tte_outcome, statuskm) ~ 1, data = testdata)
testfitstatusci <- survfit(Surv(tte_outcome, statusci) ~ 1, data = testdata)

testdatakm <- broom::tidy(testfitkm) %>%
  mutate(
    summand = n.event / ((n.risk - n.event) * n.risk),
    surv=cumprod(1 - n.event / n.risk),
    #surv.ll = conf.low,
    #surv.ul = conf.high,
    surv.se = surv * sqrt(cumsum(summand))
  )


testdatastatuskm0 <- broom::tidy(testfitstatuskm)
testdatastatuskm <-
  left_join(
    testdatastatuskm0 %>% filter(state=="(s0)") %>% select(time, n.risk, n.censor),
    testdatastatuskm0 %>% filter(state==outcome) %>% select(-n.censor, -n.risk),
    by="time"
  ) %>%
  mutate(
    inv.estimate=1-estimate, inv.conf.low=1-conf.high, inv.conf.high=1-conf.low,
    summand = n.event / ((n.risk - n.event) * n.risk),
    surv = cumprod(1 - n.event / n.risk),
    #surv.ll = conf.low,
    #surv.ul = conf.high,
    surv.se = surv * sqrt(cumsum(summand))
  )

testdatastatusci0 <- broom::tidy(testfitstatusci)
testdatastatusci <-
  bind_cols(
    testdatastatusci0 %>% filter(state=="(s0)") %>% select(time, n.risk, n.censor),
    testdatastatusci0 %>% filter(state!="(s0)") %>% group_by(time) %>% summarise(n.allevents=sum(n.event)) %>% ungroup() %>% select(-time),
    testdatastatusci0 %>% filter(state==outcome) %>% select(-time, -n.censor, -n.risk, -state),
  ) %>%
  mutate(
    risk = estimate,
    surv = 1-risk,
    risk.se = std.error,
    surv.se = std.error,
    surv.ll = conf.low,
    surv.ul = conf.high,

    kmsummand = n.allevents / ((n.risk - n.allevents) * n.risk),
    kmsurv = cumprod(1 - n.allevents / n.risk),
    summand = (n.event / n.risk) * lag(kmsurv, 1, 1),
    risk = cumsum(summand),
    kmsurv.se = kmsurv * sqrt(cumsum(kmsummand)),
    risk.se = cif.se(time, estimate, n.risk, n.event, kmsurv, kmsummand)
  )

