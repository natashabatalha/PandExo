$(document).ready(function () {
    if (!window.console)
        window.console = {};

    if (!window.console.log)
        window.console.log = function() {};

    // If we're on the main page, startup checking of the runs
    $("#calculation-table tbody tr").each(function(i, row) {
        var status_url = "/calculation/statusspec/" + row.id.slice(1);
        updateStatus(status_url);
    });
});

function updateStatus(statusUrl) {
    console.log("Updating status" + statusUrl);
    // Query the status url for information about the state of the calculation
    $.getJSON(statusUrl, function(response) {
        $('#m' + response.id).replaceWith(response.html);
        if ((response.state != 'finished') && (response.state != 'cancelled')) {
            // Check in another 2 seconds
            setTimeout(function() {
                updateStatus(statusUrl);
            }, 2000);
        }
    });
}
