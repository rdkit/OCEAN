from django import template
from ocean.views import getTC
register = template.Library()
import ocean.settings
from os import path
import os


class DataParser(template.Node):
    def __init__(self,varname=None,data=None,directShow=False):
        self.varname = varname
        self.data = data
        self.directShow = directShow

    def __repr__(self):
        return "<DataParser>"

    def render(self,context):
        if self.directShow:
            return self.data
        else:
            context[self.varname] = self.data
            return ''


@register.filter(name="chembl_link_target")
def chembl_link_target(target):
    return ocean.settings.CHEMBL_LINK_TARGET.format(target)

@register.filter(name="chembl_link_compound")
def chembl_link_compound(compound):
    return ocean.settings.CHEMBL_LINK_COMPOUND.format(compound)

@register.filter(name="tc")
def tc(thisSmiles,otherMolregno):
    return getTC(thisSmiles,otherMolregno)

@register.tag
def includeIfExist(parser,token):
    tokens = token.contents.split()
    PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
    fileCheck = "{}/../../templates/{}".format(PROJECT_ROOT,tokens[1])
    if path.isfile(fileCheck):
        with open(fileCheck,'r') as fh:
            data = fh.read()
    else:
        data = ""
    return DataParser(data=data, directShow=True)

@register.tag
def getDataSources(parser, token):
    tokens = token.contents.split()
    data = [ds.name for ds in ocean.settings.DATASOURCES]
    return DataParser(varname = tokens[2], data=data)

@register.tag
def getFingerprintIDs(parser,token):
    tokens = token.contents.split()
    data = range(len(ocean.settings.FP_METHODS))
    return DataParser(varname=tokens[2], data=data)

@register.tag
def loadAdminHelpText(parser,token):
    tokens = token.contents.split()
    data = ocean.settings.ADMIN_HELPTEXT
    return DataParser(varname=tokens[2], data=data)