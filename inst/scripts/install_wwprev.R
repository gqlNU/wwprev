
local.build <- FALSE
if (local.build) {
	#devtools::create('wwprev')
	workdir <- '/Volumes/WorkSpace/OnGitHub/wwprev/'
	setwd(workdir)
    devtools::document()
    devtools::load_all()
} else {
    PAT <- readline(prompt='Enter PAT here: ') # enter token here
    devtools::install_github("gqlNU/wwprev",auth_token = PAT,force=TRUE)
}
