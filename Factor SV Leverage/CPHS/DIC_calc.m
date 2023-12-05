%DIC calc

llh_post_used=llh_post(2001:end);
llh_store_used=llh_store(2001:end);

[~,id] = max(llh_post_used);
DIC = -4*mean(llh_store_used) + 2*llh_store_used(id);