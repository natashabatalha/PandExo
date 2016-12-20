$(document).ready(function () {
    if (!window.console)
        window.console = {};

    if (!window.console.log)
        window.console.log = function() {};

    // Attach events to the form
    $("body").on("#calculation-form", "submit", function() {
        console.log("Submitted job");
        newCalculation($(this));
        return false;
    });
});


function newCalculation(form) {
    // Submit a new calculation using AJAX
    args = form.formToDict();
    args._xsrf = getCookie("_xsrf");
    $.ajax({
        data: $.param(args),
        type: 'POST',
        url: '/calculation/newspec',
        success: function(response) {
            response = eval("(" + response + ")");
            status_url = response.location;
            $('#inbox').append($(response.html));
            updateStatus(status_url);
        },
        error: function() {
            alert('Unexpected error');
        }
    });
}

jQuery.fn.formToDict = function() {
    var fields = this.serializeArray();
    var json = {};
    for (var i = 0; i < fields.length; i++) {
        json[fields[i].name] = fields[i].value;
    }
    if (json.next) delete json.next;
    return json;
};


function getCookie(name) {
    var r = document.cookie.match("\\b" + name + "=([^;]*)\\b");
    return r ? r[1] : undefined;
}