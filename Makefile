all: notebooks/DAIR3_project.html

notebooks/DAIR3_project.html: notebooks/DAIR3_project.Rmd
	Rscript -e "rmarkdown::render('notebooks/DAIR3_project.Rmd')"

clean:
	rm -f notebooks/*.html
