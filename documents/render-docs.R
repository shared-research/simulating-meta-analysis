# render documents

# paper

rmarkdown::render(input = "documents/paper/paper.Rmd", 
                  output_dir = "documents/output")

# supplementary

rmarkdown::render(input = "documents/paper/suppl.Rmd", 
                  output_dir = "documents/output",
                  output_file = "suppl.pdf")