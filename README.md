
# Editing the Documentation:


Find the source files and further instructions in:

> master@doc/
> master@doc/README

------------------------------

# Updating the web pages:


Compile in the branch "master" the html pages:

> master@doc/: make html

Update the branch "gh-pages":

> master@.: git checkout gh-pages
> gh-pages@.: ./update_html
> gh-pages@.: git push

