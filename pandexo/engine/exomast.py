import logging
import re
import requests
import urllib

import tornado.escape
from tornado.httpclient import AsyncHTTPClient

logger = logging.getLogger(__name__)
IDENTIFIER_URL = "https://exo.mast.stsci.edu/api/v0.1/exoplanets/identifiers/"
NO_TARGET_DATA_ERROR = 'Whoops, no data for this target!'
PREFERRED_CATALOGS = ('nexsci', 'exoplanets.org')


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


def _select_preferred_catalog(target_data):
    # Some targets have multiple catalogs; nexsci is the first choice.
    if not target_data:
        raise Exception(NO_TARGET_DATA_ERROR)

    if len(target_data) > 1:
        catalog_dict = {data['catalog_name']: index for index, data in enumerate(target_data)}
        for catalog in PREFERRED_CATALOGS:
            if catalog in catalog_dict:
                return target_data[catalog_dict[catalog]]
        return target_data[0]
    return target_data[0]


def _canonical_name_from_identifier_payload(planetnames):
    return planetnames['canonicalName']


def _format_target_result(canonical_name, target_data):
    target_data = _select_preferred_catalog(target_data)
    url_name = re.sub(r'\W+', '', canonical_name)
    url = 'https://exo.mast.stsci.edu/exomast_planet.html?planet={}'.format(url_name)
    return target_data, url


def get_canonical_name(target_name):
    '''Get ExoMAST prefered name for exoplanet.
        Parameters
        ----------
        target_name : string
            The name of the target transit.
        Returns
        -------
        canonical_name : string
    '''

    params = {"name": target_name}

    response = requests.get(IDENTIFIER_URL, params=params)
    return _canonical_name_from_identifier_payload(response.json())


async def async_get_canonical_name(target_name):
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

    target_url = f"{IDENTIFIER_URL}?name={url_name}"
    logger.debug(f"\tQuerying {target_url}")

    http_client = AsyncHTTPClient()
    try:
        response = await http_client.fetch(target_url)
    except Exception as e:
        logger.error(f"\t{e}")
        raise
    else:
        canonical_name = _canonical_name_from_identifier_payload(
            tornado.escape.json_decode(response.body)
        )
        logger.debug(f"\tCanonical name is {canonical_name}")

    return canonical_name

def get_target_data(target_name):
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

    canonical_name = get_canonical_name(target_name)
    target_url = build_target_url(canonical_name)

    response = requests.get(target_url)
    if response.status_code == 200:
        target_data = response.json()
    else:
        raise Exception(NO_TARGET_DATA_ERROR)

    return _format_target_result(canonical_name, target_data)


async def async_get_target_data(target_name):
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

    canonical_name = await async_get_canonical_name(target_name)

    target_url = build_target_url(canonical_name)
    logger.debug(f"\tTarget URL is {target_url}")

    http_client = AsyncHTTPClient()
    try:
        response = await http_client.fetch(target_url)
    except Exception as e:
        logger.error(f"\t{e}")
        raise
    else:
        logger.debug(f"\tResponse is {response}")
        if response.code == 200:
            target_data = tornado.escape.json_decode(response.body)
        else:
            logger.error(f"\tResponse code was {response.code}")
            raise Exception(NO_TARGET_DATA_ERROR)

    return _format_target_result(canonical_name, target_data)
