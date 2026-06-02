import itertools
import logging
import os
import re
import sys
import urllib

import tornado.escape
from tornado.httpclient import AsyncHTTPClient

logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def build_target_url(target_name):
    '''Build restful api url based on target name.
    Parameters
        ----------
        target_name : string
            The name of the target transit.
        Returns
        -------
        target_url : string
    '''
    # Encode the target name string.
    encode_target_name = urllib.parse.quote(target_name, encoding='utf-8')
    target_url = "https://exo.mast.stsci.edu/api/v0.1/exoplanets/{}/properties/".format(encode_target_name)

    return target_url

async def get_canonical_name(target_name):
    '''Get ExoMAST prefered name for exoplanet.
        Parameters
        ----------
        target_name : string
            The name of the target transit.
        Returns
        -------
        canonical_name : string
    '''

    url_name = urllib.parse.quote_plus(target_name)

    target_url = f"https://exo.mast.stsci.edu/api/v0.1/exoplanets/identifiers/?name={url_name}"
    logging.debug(f"\tQuerying {target_url}")

    http_client = AsyncHTTPClient()
    try:
        response = await http_client.fetch(target_url)
    except Exception as e:
        logging.error(f"\t{e}")
        canonical_name = ""
    else:
        planetnames = tornado.escape.json_decode(response.body)
        canonical_name = planetnames['canonicalName']
        logging.debug(f"\tCanonical name is {canonical_name}")
    

    return canonical_name

async def get_target_data(target_name):
    """
    Send request to exomast restful api for target information.
    Parameters
    ----------
    target_name : string
        The name of the target transit
    Returns
    -------
    target_data: json:
        json object with target data.
    """

    canonical_name = await get_canonical_name(target_name)

    target_url = build_target_url(canonical_name)
    logging.debug(f"\tTarget URL is {target_url}")

    http_client = AsyncHTTPClient()
    try:
        response = await http_client.fetch(target_url)
    except Exception as e:
        logging.error(f"\t{e}")
        target_data = []
    else:
        logging.debug(f"\tResponse is {response}")
        if response.code == 200:
            target_data = tornado.escape.json_decode(response.body)
        else:
            logging.error(f"\tResponse code was {response.code}")
            raise Exception('Whoops, no data for this target!')

        # Some targets have multiple catalogs
        # nexsci is the first choice.
        if len(target_data) > 1:
            # Get catalog names from exomast and make then the keys of a dictionary
            # and the values are its position in the json object.
            catalog_dict = {data['catalog_name']: index for index, data in enumerate(target_data)}

            # Parse based on catalog accuracy.
            if 'nexsci' in list(catalog_dict.keys()):
                target_data = target_data[catalog_dict['nexsci']]
            elif 'exoplanets.org' in list(catalog_dict.keys()):
                target_data = target_data[catalog_dict['exoplanets.org']]
            else:
                target_data = target_data[0]
        else:
            target_data = target_data[0]

    # Strip spaces and non numeric or alphabetic characters and combine.
    url = 'https://exo.mast.stsci.edu/exomast_planet.html?planet={}'.format(re.sub(r'\W+', '', canonical_name))

    return target_data, url
