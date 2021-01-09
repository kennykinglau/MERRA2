import glob

currentSite = "MaunaKea"
siteInitial = "MK"
gndData = 2
cloud = "vapor"

failedJobs = []

for i in glob.glob("/n/home10/anwang16/merra2_work/merra2_products/%s/farmfiles/jobs/*"%(currentSite)):
    picklefile = i.replace(".sh","_gndData2.p")
    failedJobs.append(picklefile)

print failedJobs
