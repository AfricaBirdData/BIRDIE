
# In case we need to install dependencies we can do something like this:
# withr::with_libpaths(new="/home/crvfra001/R/x86_64-pc-linux-gnu-library/4.2",
#                      install.packages("rnaturalearth", repos = "https://cran.mirror.ac.za/"))

withr::with_libpaths(new="/home/crvfra001/R/x86_64-pc-linux-gnu-library/4.2",
                     remotes::install_github('AfricaBirdData/BIRDIE',
                                             auth_token = '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',
                                             dependencies = FALSE) )

withr::with_libpaths(new="/home/crvfra001/R/x86_64-pc-linux-gnu-library/4.2",
                     remotes::install_github('AfricaBirdData/ABAP',
                                             auth_token = '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',
                                             dependencies = FALSE) )

withr::with_libpaths(new="/home/crvfra001/R/x86_64-pc-linux-gnu-library/4.2",
                     remotes::install_github('AfricaBirdData/CWAC',
                                             auth_token = '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',
                                             dependencies = FALSE) )

withr::with_libpaths(new="/home/crvfra001/R/x86_64-pc-linux-gnu-library/4.2",
                     remotes::install_github('AfricaBirdData/ABDtools',
                                             auth_token = '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',
                                             dependencies = FALSE) )
