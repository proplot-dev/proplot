/* Toggle light and dark theme */
/* See: https://dev.to/ananyaneogi/create-a-dark-light-mode-switch-with-css-variables-34l8 */
/* CSS files from: https://github.com/richleland/pygments-css */
/* Toggle functions and event handler */
const toggleSwitch = document.getElementById('lightdark-checkbox');
const pygmentsLink = document.getElementById('pygments-style');
function lightToggle() {
    document.documentElement.setAttribute('data-theme', 'light');
    pygmentsLink.href = 'https://rawgit.com/richleland/pygments-css/master/friendly.css';
    localStorage.setItem('theme', 'light');
}
function darkToggle() {
    document.documentElement.setAttribute('data-theme', 'dark');
    pygmentsLink.href = 'https://rawgit.com/richleland/pygments-css/master/monokai.css';
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
