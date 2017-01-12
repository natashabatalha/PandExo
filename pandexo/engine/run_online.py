import json
import os
import uuid
from collections import namedtuple, OrderedDict
import tornado.escape
import tornado.httpserver
import tornado.ioloop
import tornado.options
import tornado.web
from concurrent.futures import ProcessPoolExecutor
from pandexo import wrapper
from tornado.options import define, options
import pickle
from ComputeZ import computeAlpha
from utils.plotters import create_component_jwst, create_component_spec, create_component_hst
import pandas as pd 
import numpy as np

#define location of temp files
__TEMP__ = os.path.join(os.path.dirname(__file__), "temp")

define("port", default=1111, help="run on the given port", type=int)

# Define a simple named tuple to keep track for submitted calculations
CalculationTask = namedtuple('CalculationTask', ['id', 'name', 'task',
                                                 'cookie', 'count'])


class Application(tornado.web.Application):
    """Gobal settings of the server
    This defines the global settings of the server. This parses out the
    handlers, and includes settings for if ever we want to tie this to a
    database.
    """
    def __init__(self):
        handlers = [
            (r"/", HomeHandler),
            (r"/about", AboutHandler),
            (r"/dashboard", DashboardHandler),
            (r"/dashboardspec", DashboardSpecHandler),
            (r"/dashboardhst", DashboardHSTHandler),
            (r"/tables", TablesHandler),
            (r"/helpfulplots", HelpfulPlotsHandler),
            (r"/calculation/new", CalculationNewHandler),
            (r"/calculation/newHST", CalculationNewHSTHandler),
            (r"/calculation/newspec", CalculationNewSpecHandler),
            (r"/calculation/status/([^/]+)", CalculationStatusHandler),
            (r"/calculation/statushst/([^/]+)", CalculationStatusHSTHandler),
            (r"/calculation/statusspec/([^/]+)", CalculationStatusSpecHandler),
            (r"/calculation/view/([^/]+)", CalculationViewHandler),
            (r"/calculation/viewhst/([^/]+)", CalculationViewHSTHandler),
            (r"/calculation/viewspec/([^/]+)", CalculationViewSpecHandler),
            (r"/calculation/download/([^/]+)", CalculationDownloadHandler),
            (r"/calculation/downloadspec/([^/]+)", CalculationDownloadSpecHandler),
            (r"/calculation/downloadpandin/([^/]+)", CalculationDownloadPandInHandler)
        ]
        settings = dict(
            blog_title=u"Pandexo",
            template_path=os.path.join(os.path.dirname(__file__), "templates"),
            static_path=os.path.join(os.path.dirname(__file__), "static"),
            xsrf_cookies=True,
            cookie_secret="__TODO:_GENERATE_YOUR_OWN_RANDOM_VALUE_HERE__",
            debug=False,
        )
        super(Application, self).__init__(handlers, **settings)


class BaseHandler(tornado.web.RequestHandler):
    """
    Logic to handle user information and database access might go here.
    """
    executor = ProcessPoolExecutor(max_workers=4)
    buffer = OrderedDict()

    def _get_task_response(self, id):
        """
        Simple function to grab a calculation that's stored in the buffer,
        and return a dictionary/json-like response to the front-end.
        """
        calc_task = self.buffer.get(id)
        task = calc_task.task

        response = {'id': id,
                    'name': calc_task.name,
                    'count': calc_task.count}

        if task.running():
            response['state'] = 'running'
            response['code'] = 202
        elif task.done():
            response['state'] = 'finished'
        elif task.cancelled():
            response['state'] = 'cancelled'
        else:
            response['state'] = 'pending'

        response['html'] = tornado.escape.to_basestring(
            self.render_string("calc_row.html", response=response))
        return response

    def _get_task_response_spec(self, id):
        """
        Simple function to grab a calculation that's stored in the buffer,
        and return a dictionary/json-like response to the front-end.
        """
        calc_task = self.buffer.get(id)
        task = calc_task.task

        response = {'id': id,
                    'name': calc_task.name,
                    'count': calc_task.count}

        if task.running():
            response['state'] = 'running'
            response['code'] = 202
        elif task.done():
            response['state'] = 'finished'
        elif task.cancelled():
            response['state'] = 'cancelled'
        else:
            response['state'] = 'pending'

        response['html'] = tornado.escape.to_basestring(
            self.render_string("calc_rowspec.html", response=response))

        return response

    def _get_task_response_hst(self, id):
        """
        Simple function to grab a calculation that's stored in the buffer,
        and return a dictionary/json-like response to the front-end.
        """
        calc_task = self.buffer.get(id)
        task = calc_task.task

        response = {'id': id,
                    'name': calc_task.name,
                    'count': calc_task.count}

        if task.running():
            response['state'] = 'running'
            response['code'] = 202
        elif task.done():
            response['state'] = 'finished'
        elif task.cancelled():
            response['state'] = 'cancelled'
        else:
            response['state'] = 'pending'

        response['html'] = tornado.escape.to_basestring(
            self.render_string("calc_rowhst.html", response=response))

        return response
        
    def write_error(self, status_code, **kwargs):
        """
        This renders a customized error page
        """
        self.render('errors.html',page=None)



    def _get_task_result(self, id):
        """
        This method grabs only the result returned from the python `Future`
        object. This contains the stuff that Pandeia returns.
        """
        calc_task = self.buffer.get(id)
        task = calc_task.task

        return task.result()

    def _add_task(self, id, name, task):
        """
        This creates the task and adds it to the buffer.
        """
        self.buffer[id] = CalculationTask(id=id, name=name, task=task,
                                          count=len(self.buffer)+1,
                                          cookie=self.get_cookie("pandexo_user"))

        # Only allow 100 tasks **globally**. This will delete old tasks first.
        if len(self.buffer) > 15:
            self.buffer.popitem(last=False)
            


class HomeHandler(BaseHandler):
    def get(self):
        """
        This sets an **unsecured** cookie. If user accounts gets
        implemented, this must be changed to a secure cookie.
        """
        if not self.get_cookie("pandexo_user"):
            self.set_cookie("pandexo_user", str(uuid.uuid4()))

        self.render("home.html")
        
class AboutHandler(BaseHandler):
    def get(self):
        """
        Render about PandExo Page
        """
        self.render("about.html")

class TablesHandler(BaseHandler):
    def get(self):
        """
        Render tables with confirmed candidates
        """
        self.render("tables.html")

class HelpfulPlotsHandler(BaseHandler):
    def get(self):
        """
        Renders helpful bokeh plots
        """
        self.render("helpfulplots.html")


class DashboardHandler(BaseHandler):
    """
    Request handler for the dashboard page. This will retrieve and render
    the html template, along with the list of current task objects.
    """
    def get(self):
        task_responses = [self._get_task_response(id) for id, nt in
                          self.buffer.items()
                          if ((nt.cookie == self.get_cookie("pandexo_user"))
                          & (id[len(id)-1]=='e'))]
        
        self.render("dashboard.html", calculations=task_responses[::-1])

class DashboardHSTHandler(BaseHandler):
    """
    Request handler for the dashboard page. This will retrieve and render
    the html template, along with the list of current task objects.
    """
    def get(self):
        task_responses = [self._get_task_response_hst(id) for id, nt in
                          self.buffer.items()
                          if ((nt.cookie == self.get_cookie("pandexo_user"))
                          & (id[len(id)-1]=='h'))]
        
        self.render("dashboardhst.html", calculations=task_responses[::-1])



class DashboardSpecHandler(BaseHandler):
    """
    Request handler for the dashboard page. This will retrieve and render
    the html template, along with the list of current task objects.
    """
    def get(self):
        task_responses = [self._get_task_response_spec(id) for id, nt in
                          self.buffer.items()
                          if ((nt.cookie == self.get_cookie("pandexo_user"))
                          & (id[len(id)-1]=='s'))]
        self.render("dashboardspec.html", calculations=task_responses[::-1])


class CalculationNewHandler(BaseHandler):
    """
    This request handler deals with processing the form data and submitting
    a new calculation task to the parallelized workers.
    """
    def get(self):
        self.render("new.html", id=id)

    def post(self):
        """
        The post method contains the retured data from the form data (
        accessed by using `self.get_argument(...)` for specific arguments,
        or `self.request.body` to grab the entire returned object.
        """
        
        #print(self.request.body)
        
        id = str(uuid.uuid4())+'e'
                
        #upload planet file
        fileinfo_plan = self.request.files['planFile'][0]
        fname_plan = fileinfo_plan['filename']
        extn_plan = os.path.splitext(fname_plan)[1]
        cname_plan = id+'planet' + extn_plan
        fh_plan = open(os.path.join(__TEMP__, cname_plan), 'w')
        fh_plan.write(fileinfo_plan['body'])
        
        with open(os.path.join(os.path.dirname(__file__), "reference",
                               "exo_input.json")) as data_file:

            exodata = json.load(data_file)
            exodata["telescope"] = 'jwst'
            exodata["calculation"] = 'fml' #always for online form
            exodata["star"]["type"] = self.get_argument("type")
            if exodata["star"]["type"] == "user":     
                fileinfo_star = self.request.files['starFile'][0]
                fname_star = fileinfo_star['filename']
                extn_star = os.path.splitext(fname_star)[1]
                cname_star = id+'star' + extn_star
                fh_star = open(os.path.join(__TEMP__, cname_star), 'w')
                fh_star.write(fileinfo_star['body'])
                exodata["star"]["starpath"] = os.path.join(__TEMP__, cname_star)
                exodata["star"]["f_unit"] = self.get_argument("starfunits")
                exodata["star"]["w_unit"] = self.get_argument("starwunits")
            else: 
                exodata["star"]["temp"] = float(self.get_argument("temp"))
                exodata["star"]["logg"] = float(self.get_argument("logg"))
                exodata["star"]["metal"] = float(self.get_argument("metal"))
            exodata["star"]["mag"] = float(self.get_argument("mag"))
            exodata["star"]["ref_wave"] = float(self.get_argument("ref_wave"))
            exodata["planet"]["exopath"] = os.path.join(__TEMP__, cname_plan)
            exodata["planet"]["w_unit"] = self.get_argument("planwunits")
            exodata["planet"]["f_unit"] = self.get_argument("planfunits")
            exodata["observation"]["fraction"] = float(self.get_argument("fraction"))
            exodata["observation"]["noccultations"] = float(self.get_argument("numtrans"))
            exodata["observation"]["sat_level"] = float(self.get_argument("satlevel"))
            
            #for phase curves user doen't necessarily have to input a transit duration 
            try:
                exodata["planet"]["transit_duration"] = float(self.get_argument("transit_duration"))
            except:
                #but if they dont.. make sure that the planet units are in seconds... 
                if exodata["planet"]["w_unit"] == 'sec':
                    exodata["planet"]["transit_duration"] = 0.0
                else: 
                    print "Need to give transit duration"
                    raise 
                
            #noise floor, set to 0.0 of no values are input        
            try: 
                exodata["observation"]["noise_floor"] = float(self.get_argument("noisefloor"))
            except:
                exodata["observation"]["noise_floor"] = 0.0
            try: 
                fileinfo_noise = self.request.files['noisefile'][0]
                fname_noise = fileinfo_noise['filename']
                extn_noise = os.path.splitext(fname_noise)[1]
                cname_noise = id+'noise' + extn_noise
                fh_noise = open(os.path.join(__TEMP__, cname_noise), 'w')
                fh_noise.write(fileinfo_star['body'])
                exodata["observation"]["noise_floor"] = os.path.join(__TEMP__, cname_noise)
            except:
                exodata["observation"]["noise_floor"] = 0.0
                
        if (self.get_argument("instrument")=="MIRI"): 
            with open(os.path.join(os.path.dirname(__file__), "reference",
                               "miri_input.json")) as data_file:   
                pandata = json.load(data_file)       
                mirimode = self.get_argument("mirimode")
                if (mirimode == "lrsslit"):
                    pandata["configuration"]["mode"] = mirimode
                    pandata["configuration"]["instrument"]["aperture"]="lrsslit"
                    
        if (self.get_argument("instrument")=="NIRSpec"): 
            with open(os.path.join(os.path.dirname(__file__), "reference",
                               "nirspec_input.json")) as data_file:
                pandata = json.load(data_file)  
                nirspecmode = self.get_argument("nirspecmode")
                pandata["configuration"]["instrument"]["disperser"] = nirspecmode[0:5]
                pandata["configuration"]["instrument"]["filter"] = nirspecmode[5:11]
                pandata["configuration"]["detector"]["subarray"]  = self.get_argument("nirspecsubarray")
        if (self.get_argument("instrument")=="NIRCam"): 
            with open(os.path.join(os.path.dirname(__file__), "reference",
                               "nircam_input.json")) as data_file:
                pandata = json.load(data_file) 
                pandata["configuration"]["instrument"]["filter"]  = self.get_argument("nircammode")
                pandata["configuration"]["detector"]["subarray"]  = self.get_argument("nircamsubarray")
        if (self.get_argument("instrument")=="NIRISS"): 
            with open(os.path.join(os.path.dirname(__file__), "reference",
                               "niriss_input.json")) as data_file:
                pandata = json.load(data_file) 
                nirissmode = self.get_argument("nirissmode")
                pandata["configuration"]["detector"]["subarray"] = nirissmode
        
        #write in optimal groups or set a number 
        try:
            pandata["configuration"]["detector"]["ngroup"] = int(self.get_argument("optimize"))
        except: 
            pandata["configuration"]["detector"]["ngroup"] = self.get_argument("optimize")
        
        finaldata = {"pandeia_input": pandata , "pandexo_input":exodata}

        task = self.executor.submit(wrapper, finaldata)


        self._add_task(id, self.get_argument("calcName"), task)

        response = self._get_task_response(id)
        response['info'] = {}
        response['location'] = '/calculation/status/{}'.format(id)
        
        
        self.write(dict(response))
        self.redirect("/dashboard")
        
    
class CalculationNewHSTHandler(BaseHandler):
    """
    This request handler deals with processing the form data and submitting
    a new HST calculation task to the parallelized workers.
    """
    def get(self):
        self.render("newHST.html", id=id)

    def post(self):
        """
        The post method contains the retured data from the form data (
        accessed by using `self.get_argument(...)` for specific arguments,
        or `self.request.body` to grab the entire returned object.
        """
        
        #print(self.request.body)
        
        id = str(uuid.uuid4())+'h'
                
        #upload planet file
        fileinfo_plan = self.request.files['planFile'][0]
        fname_plan = fileinfo_plan['filename']
        extn_plan = os.path.splitext(fname_plan)[1]
        cname_plan = id+'planet' + extn_plan
        fh_plan = open(os.path.join(__TEMP__, cname_plan), 'w')
        fh_plan.write(fileinfo_plan['body'])
        
        with open(os.path.join(os.path.dirname(__file__), "reference",
                               "exo_input.json")) as data_file:
            exodata = json.load(data_file)
            exodata["telescope"] = 'hst'
            exodata["star"]["mag"]         = float(self.get_argument("mag"))
            exodata["star"]["ref_wave"]     = float(self.get_argument("ref_wave"))
            exodata["planet"]["exopath"]    = os.path.join(__TEMP__, cname_plan)
            exodata["planet"]["w_unit"]     = self.get_argument("planwunits")
            exodata["planet"]["f_unit"]     = self.get_argument("planfunits")
            exodata["planet"]["depth"]      = float(self.get_argument("depth"))
            exodata["planet"]["i"]          = float(self.get_argument("i"))
            exodata["planet"]["ars"]        = float(self.get_argument("ars"))
            exodata["planet"]["period"]     = float(self.get_argument("period"))
            exodata["planet"]["ecc"]        = float(self.get_argument("ecc"))
            try:
                exodata["planet"]["w"]      = float(self.get_argument("w"))
            except:
                exodata["planet"]["w"]      = 90.
            exodata["planet"]["transit_duration"]   = float(self.get_argument("transit_duration"))
            exodata["observation"]["noise_floor"]           = 0.0
            exodata["calculation"]                          = 'scale'
                
        if (self.get_argument("instrument")=="STIS"): 
            with open(os.path.join(os.path.dirname(__file__), "reference",
                               "stis_input.json")) as data_file:   
                pandata = json.load(data_file)       
                stismode = self.get_argument("stismode")
        if (self.get_argument("instrument")=="WFC3"): 
            with open(os.path.join(os.path.dirname(__file__), "reference",
                               "wfc3_input.json")) as data_file:
                pandata = json.load(data_file)  
                pandata["configuration"]['detector']['subarray']    = self.get_argument("subarray")
                pandata["configuration"]['detector']['nsamp']       = int(self.get_argument("nsamp"))
                pandata["configuration"]['detector']['samp_seq']    = self.get_argument("samp_seq")
                pandata["configuration"]['instrument']['disperser'] = self.get_argument("wfc3mode")
            try: 
                pandata["strategy"]["norbits"]           = int(self.get_argument("norbits"))
            except:
                pandata["strategy"]["norbits"]           = None
            exodata["observation"]["noccultations"]         = int(self.get_argument("noccultations"))
            pandata["strategy"]["nchan"]                 = int(self.get_argument("nchan"))
            pandata["strategy"]["scanDirection"]         = self.get_argument("scanDirection")
            try:
                pandata["strategy"]["windowSize"]        = float(self.get_argument("windowSize"))
            except:
                pandata["strategy"]["windowSize"]        = 20.
            pandata["strategy"]["schedulability"]           = self.get_argument("schedulability")

        finaldata = {"pandeia_input": pandata , "pandexo_input":exodata}

        task = self.executor.submit(wrapper, finaldata)

        self._add_task(id, self.get_argument("calcName"), task)

        response = self._get_task_response_hst(id)
        response['info'] = {}
        response['location'] = '/calculation/statushst/{}'.format(id)
        
        
        self.write(dict(response))
        self.redirect("/dashboardhst")
        
            
class CalculationNewSpecHandler(BaseHandler):
    """
    This request handler deals with processing the form data and submitting
    a new calculation task to the parallelized workers.
    """
    def get(self):
        self.render("newspec.html", id=id)

    def post(self):
        """
        The post method contains the retured data from the form data (
        accessed by using `self.get_argument(...)` for specific arguments,
        or `self.request.body` to grab the entire returned object.
        """
        
        #print(self.request.body)
        
        id = str(uuid.uuid4())+'s'
        finaldata= {}      
        mols = {}
        mols['H2'] = 0.0
        mols['He'] = 0.0
        try: 
            mols[self.get_argument("mol1")] = float(self.get_argument("mol1_vmr"))
        except:
            if self.get_argument("mol1_vmr").replace(' ','').lower() == 'bkg': 
                mols[self.get_argument("mol1")] = self.get_argument("mol1_vmr")
        try: 
            mols[self.get_argument("mol2")] = float(self.get_argument("mol2_vmr"))
        except:
            if self.get_argument("mol2_vmr").replace(' ','').lower() == 'bkg': 
                mols[self.get_argument("mol2")] = self.get_argument("mol2_vmr")
        try: 
            mols[self.get_argument("mol3")] = float(self.get_argument("mol3_vmr"))
        except:
            if self.get_argument("mol3_vmr").replace(' ','').lower() == 'bkg': 
                mols[self.get_argument("mol3")] = self.get_argument("mol3_vmr")
        try: 
            mols[self.get_argument("mol4")] = float(self.get_argument("mol4_vmr"))
        except:
            if self.get_argument("mol4_vmr").replace(' ','').lower() == 'bkg': 
                mols[self.get_argument("mol4")] = self.get_argument("mol4_vmr")
        try: 
            mols[self.get_argument("mol5")] = float(self.get_argument("mol5_vmr"))
        except:
            if self.get_argument("mol5_vmr").replace(' ','').lower() == 'bkg': 
                mols[self.get_argument("mol5")] = self.get_argument("mol5_vmr")
        try: 
            mols[self.get_argument("mol6")] = float(self.get_argument("mol6_vmr"))
        except:
            if self.get_argument("mol6_vmr").replace(' ','').lower() == 'bkg': 
                mols[self.get_argument("mol6")] = self.get_argument("mol6_vmr")
        try: 
            mols[self.get_argument("mol7")] = float(self.get_argument("mol7_vmr"))
        except:
            if self.get_argument("mol7_vmr").replace(' ','').lower() == 'bkg': 
                mols[self.get_argument("mol7")] = self.get_argument("mol7_vmr")
        try: 
            mols[self.get_argument("mol8")] = float(self.get_argument("mol8_vmr"))
        except:
            if self.get_argument("mol8_vmr").replace(' ','').lower() == 'bkg': 
                mols[self.get_argument("mol8")] = self.get_argument("mol8_vmr")
        try: 
            mols[self.get_argument("mol9")] = float(self.get_argument("mol9_vmr"))
        except:
            if self.get_argument("mol9_vmr").replace(' ','').lower() == 'bkg': 
                mols[self.get_argument("mol9")] = self.get_argument("mol9_vmr")
        try: 
            mols[self.get_argument("mol10")] = float(self.get_argument("mol10_vmr"))
        except:
            if self.get_argument("mol10_vmr").replace(' ','').lower() == 'bkg': 
                mols[self.get_argument("mol10")] = self.get_argument("mol10_vmr")
        

        finaldata["mols"] = mols
        finaldata['T'] = float(self.get_argument("T"))
        finaldata['g'] = float(self.get_argument("g"))
        finaldata['Rp'] = float(self.get_argument("Rp"))
        finaldata['R*'] = float(self.get_argument("R*"))
        finaldata['P'] = float(self.get_argument("P"))
        finaldata['P0'] = float(self.get_argument("P0"))
        try: 
            finaldata['fracH2He'] = float(self.get_argument("fracH2He"))
        except: 
            finaldata['fracH2He'] = 0.0
        finaldata['taueq'] = 0.56
        task = self.executor.submit(computeAlpha, finaldata)

        self._add_task(id, self.get_argument("calcName"), task)

        response = self._get_task_response_spec(id)
        response['info'] = {}
        response['location'] = '/calculation/statusspec/{}'.format(id)
        
        
        self.write(dict(response))
        self.redirect("/dashboardspec")


class CalculationStatusHandler(BaseHandler):
    """
    Handlers returning the status of a particular JWST calculation task.
    """
    def get(self, id):
        response = self._get_task_response(id)

        if self.request.connection.stream.closed():
            return

        self.write(dict(response))

class CalculationStatusSpecHandler(BaseHandler):
    """
    Handlers returning the status of a particular Spec calculation task.
    """
    def get(self, id):
        response = self._get_task_response_spec(id)

        if self.request.connection.stream.closed():
            return

        self.write(dict(response)) 
        
class CalculationStatusHSTHandler(BaseHandler):
    """
    Handlers returning the status of a particular HST calculation task.
    """
    def get(self, id):
        response = self._get_task_response_hst(id)

        if self.request.connection.stream.closed():
            return

        self.write(dict(response))                

class CalculationDownloadHandler(BaseHandler):
    """
    Handlers returning the downloaded data of a particular calculation task.
    Handlers returning the status of a particular calculation task.
    """
    def get(self, id):
        result = self._get_task_result(id)
  
        if self.request.connection.stream.closed():
            return
        file_name = "ETC-calculation" +id+".p"
 
        with open(os.path.join(__TEMP__,file_name), "w") as f:
            pickle.dump(result, f)
 
        buf_size = 4096
        self.set_header('Content-Type', 'application/octet-stream')
        self.set_header('Content-Disposition',
                        'attachment; filename=' + file_name)
 
        with open(os.path.join(__TEMP__,file_name), "rb") as f:
            while True:
                data = f.read(buf_size)
                if not data:
                    break
                self.write(data)
        
        allfiles = os.listdir(__TEMP__)
        for i in allfiles:
            if i.find(id) != -1:
                os.remove(os.path.join(__TEMP__,i))
        self.finish()

class CalculationDownloadSpecHandler(BaseHandler):
    """
    Handlers returning the downloaded data of a particular calculation task.
    Handlers returning the status of a particular calculation task.
    """
    def get(self, id):
        result = self._get_task_result(id)
  
        if self.request.connection.stream.closed():
            return
        file_name = "spec-calculation" +id+".p"
 
        with open(os.path.join(__TEMP__,file_name), "w") as f:
            pickle.dump(result, f)
 
        buf_size = 4096
        self.set_header('Content-Type', 'application/octet-stream')
        self.set_header('Content-Disposition',
                        'attachment; filename=' + file_name)
 
        with open(os.path.join(__TEMP__,file_name), "rb") as f:
            while True:
                data = f.read(buf_size)
                if not data:
                    break
                self.write(data)
        
        allfiles = os.listdir(__TEMP__)
        for i in allfiles:
            if i.find(id) != -1:
                os.remove(os.path.join(__TEMP__,i))
        self.finish()

class CalculationDownloadPandInHandler(BaseHandler):
    """
    Handlers returning the downloaded data of a particular calculation task.
    Handlers returning the status of a particular calculation task.
    """
    def get(self, id):
        result = self._get_task_result(id)
  
        if self.request.connection.stream.closed():
            return
        file_name = "PandExo-Input-file"+id+".txt"
 
        #with open(os.path.join(__TEMP__,file_name), "w") as f:
        #    pickle.dump(result, f)
        np.savetxt(os.path.join(__TEMP__,file_name), np.transpose([result['w'], result['alpha']]))
        
        buf_size = 4096
        self.set_header('Content-Type', 'application/octet-stream')
        self.set_header('Content-Disposition',
                        'attachment; filename=' + file_name)
 
        with open(os.path.join(__TEMP__,file_name), "rb") as f:
            while True:
                data = f.read(buf_size)
                if not data:
                    break
                self.write(data)
        
        allfiles = os.listdir(__TEMP__)
        for i in allfiles:
            if i.find(id) != -1:
                os.remove(os.path.join(__TEMP__,i))


        self.finish()



class CalculationViewHandler(BaseHandler):
    """
    This handler deals with passing the results from Pandeia to the
    `create_component_jwst` function which generates the Bokeh interative plots.
    """
    def get(self, id):
        
        result = self._get_task_result(id)
        
        script, div = create_component_jwst(result)
        div['timing_div'] = result['timing_div']
        div['input_div'] = result['input_div'] 
        div['warnings_div'] = result['warnings_div']

        #delete files
        allfiles = os.listdir(__TEMP__)
        for i in allfiles:
            if i.find(id) != -1:
                os.remove(os.path.join(__TEMP__,i))

        self.render("view.html", script=script, div=div, id=id)

class CalculationViewSpecHandler(BaseHandler):
    """
    This handler deals with passing the results from Pandeia to the
    `create_component_spec` function which generates the Bokeh interative plots.
    """
    def get(self, id):
        result = self._get_task_result(id)
        script, div = create_component_spec(result)

        self.render("viewspec.html", script=script, div=div, id=id)

class CalculationViewHSTHandler(BaseHandler):
    """
    This handler deals with passing the results from Pandeia to the
    `create_component_hst` function which generates the Bokeh interative plots.
    """
    def get(self, id):
        result = self._get_task_result(id)
        script, div = create_component_hst(result)
        div['info_div'] = result['info_div']
        self.render("viewhst.html", script=script, div=div, id=id)


def main():
    tornado.options.parse_command_line()
    http_server = tornado.httpserver.HTTPServer(Application())
    http_server.listen(options.port)
    tornado.ioloop.IOLoop.current().start()


if __name__ == "__main__":
    main()
