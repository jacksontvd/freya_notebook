import numpy as np
import matplotlib.pyplot as plt


################################################################################
# readJSONDataFile
################################################################################
import json

#-- Read JSON-formatted data file
def readJSONDataFile (filename):
  
  with open(filename) as jsonFile:
    jsonData = json.load(jsonFile)
    
    exp = list()
    c=0
    for item in jsonData.get("entries"):
      exp.append(item)
      c=c+1
    
  return exp

################################################################################
# plotExperimentalData
################################################################################

#-- Plot experimental data sets corresponding to a particular quantity
def plotExperimentalData (experiments, quantity, **keyword_parameters):
  if ('format' in keyword_parameters):
    fmtplot = keyword_parameters['format']
  else:
    fmtplot = 'ko--'
  
  checkAuthor = False
  if ('author' in keyword_parameters):
    author = keyword_parameters['author']
    checkAuthor = True
  for exp in experiments:
    if (exp['quantity']==quantity):
      if (checkAuthor):
        if (author not in exp['authors']):
          continue
      data = np.asarray(exp['data'])
      if ('renorm' in exp):
        coef = exp['renorm']
      else:
        coef = 1.0
      x = data[:,0]
      y = data[:,1]
      plt.plot(x,y*coef,fmtplot,alpha=0.5,label=exp['label'])
  return


################################################################################
# listExperimentalData
################################################################################
#-- Lists all experiments read in
def listExperimentalData (experiments):
  for exp in experiments:
    print ("{0:15s} |  {1}, {2}".format(exp['quantity'], exp['authors'], exp['year']))



################################################################################
# listExperimentalData
################################################################################
class ListTable(list):
  """ Overridden list class which takes a 2-dimensional list of
    the form [[1,2,3],[4,5,6]], and renders an HTML Table in
    IPython Notebook. """
      
  def _repr_html_(self):
    html = ["<table width=60%>"]
    html.append("<center>")
    for row in self:
      html.append("<tr>")
      for col in row:
        html.append("<td>{0}</td>".format(col))
      html.append("</tr>")
    html.append("</center>")
    return ''.join(html)

