fit1 <- mr_ash(X,y,s,s0,w0,b,20,update.w0 = FALSE)
fit2 <- mr_ash(X,y,s,s0,w0,b,20,update.w0 = TRUE)
print(max(fit1$elbo),digits = 16)
print(max(fit2$elbo),digits = 16)
