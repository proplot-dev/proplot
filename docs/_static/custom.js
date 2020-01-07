/* Toggle light and dark theme
 * See: https://dev.to/ananyaneogi/create-a-dark-light-mode-switch-with-css-variables-34l8
 * CSS files from: https://github.com/richleland/pygments-css
 * Toggle functions and event handler
 * Best light themes: pastie, friendly, murhpy
 * Best dark themes: monokai, paraiso-dark
 */
var regex = /(.*)\/.*(\.css$)/i
const toggleSwitch = document.getElementById('lightdark-checkbox');
const pygmentsLink = document.getElementById('pygments-style');
function lightToggle() {
    document.documentElement.setAttribute('data-theme', 'light');
    pygmentsLink.href = pygmentsLink.href.replace(regex, '$1/pastie$2');
    localStorage.setItem('theme', 'light');
}
function darkToggle() {
    document.documentElement.setAttribute('data-theme', 'dark');
    pygmentsLink.href = pygmentsLink.href.replace(regex, '$1/monokai$2');
    localStorage.setItem('theme', 'dark');
}
function switchTheme(e) {
    e.target.checked ? darkToggle() : lightToggle()
}
toggleSwitch.addEventListener('change', switchTheme, null);

/* Check for user preference on load */
const currentTheme = localStorage.getItem('theme') || 'light';
if (currentTheme === 'dark') {
    darkToggle();
    toggleSwitch.checked = true;
}
else {
    lightToggle();
    toggleSwitch.checked = false;
}
