const checkbox = document.getElementById('toggleCheckbox');
const content = document.getElementById('hiddenContent');

checkbox.addEventListener('change', () => 
{
    content.style.display = checkbox.checked ? 'block' : 'none';
});