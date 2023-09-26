# render documents

# paper

rmarkdown::render(input = "documents/paper/paper.Rmd", 
                  output_dir = "documents/output")

# supplementary

rmarkdown::render(input = "documents/paper/suppl.Rmd", 
                  output_dir = "documents/output",
                  output_file = "suppl.pdf")

# preprint

create_preprint()

# submission

# ## copy the paper
# fs::file_copy("documents/paper/paper.Rmd",
#               "documents/submission/paper.Rmd")
# 
# ## copy the cache
# fs::dir_copy("documents/paper/paper_cache/",
#              "documents/submission/")

rmarkdown::render(input = "documents/submission/paper.Rmd",
                  params = list(submission = TRUE))